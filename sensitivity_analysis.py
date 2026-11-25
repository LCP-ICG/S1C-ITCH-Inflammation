"""
Standalone sensitivity analysis for the 23-variable continuous-time Boolean model
of TNF-NF-kB-A20-ITCH inflammatory regulation.

Runs with Python 3 standard library only (zero external dependencies required).
Reproduces all 86 conditions x 5,000 trajectories reported in Table 3 and Table S7:
  A. Sensitivity to A20 complex logical rules (strict AND, RNF11 dispensable,
     TAX1BP1/RNF11 OR, partial ITCH-independent basal A20 activity, ITCH dispensable).
  B. Sensitivity to IkBa recovery modes (Variant D, Variant R, intermediate stochastic
     resynthesis under active IKK, Variant M).
  C. Sensitivity to relative transition rates (uniform baseline, fast signaling 5x,
     slow transcription/turnover 0.2x, asymmetric deactivation 2x / 0.5x).
  D. Sensitivity to ITCH switching frequency at fixed stationary active fraction pi.

Usage:
    python3 sensitivity_analysis.py
"""
import json
import math
import random
import time
from pathlib import Path

NODES = [
    'STIMULUS', 'ITCH_supply', 'TAX1BP1_available', 'RNF11_available', 'NFkB_import',  # 0..4 fixed
    'ITCH_functional', 'A20_mRNA', 'A20_protein', 'A20_complex',                        # 5..8
    'TNFR_drive', 'TNFR1_active', 'TRADD_recruited', 'RIPK1_signal', 'TAK1_active',     # 9..13
    'IKK_active', 'IkB_recovery', 'IkB_available',                                      # 14..16
    'NFkB_released', 'NFkB_nuclear', 'NFkB_active',                                     # 17..19
    'NFKBIA_mRNA', 'TNF_mRNA', 'TNF_secreted',                                          # 20..22
]


def wilson(k: int, n: int, z: float = 1.959963984540054) -> list[float]:
    p = k / n
    denom = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / denom
    half = z * math.sqrt((p * (1.0 - p) + z * z / (4.0 * n)) / n) / denom
    return [max(0.0, center - half), min(1.0, center + half)]


def eval_targets(s, a20_rule: str = 'strict_and', ikb_mode: str = 'D'):
    (stim, itch_sup, tax, rnf, imp,
     itch, a20_m, a20_p, a20_c,
     tnfr_d, tnfr1, tradd, ripk1, tak1,
     ikk, ikb_rec, ikb,
     nf_rel, nf_nuc, nf_act,
     ikb_m, tnf_m, tnf_s) = s

    if a20_rule == 'strict_and':
        f_a20_c = 1 if (a20_p and itch and tax and rnf) else 0
    elif a20_rule == 'no_rnf11':
        f_a20_c = 1 if (a20_p and itch and tax) else 0
    elif a20_rule == 'scaffold_or':
        f_a20_c = 1 if (a20_p and itch and (tax or rnf)) else 0
    elif a20_rule == 'itch_or_scaffolds':
        f_a20_c = 1 if (a20_p and (itch or (tax and rnf))) else 0
    else:
        raise ValueError(f"Unknown a20_rule: {a20_rule}")

    f_ikb_rec = 0 if ikk else 1
    if ikb_mode == 'D':
        f_ikb = 1 if ((ikb_rec or ikb_m) and not ikk) else 0
    elif ikb_mode == 'R':
        f_ikb = 1 if (ikb_rec or ikb_m) else 0
    elif ikb_mode == 'M_induced_escapes':
        f_ikb = 1 if (ikb_m or (ikb_rec and not ikk)) else 0
    else:
        raise ValueError(f"Unknown ikb_mode: {ikb_mode}")

    return (
        stim, itch_sup, tax, rnf, imp,
        itch_sup,                                  # ITCH_functional
        nf_act,                                    # A20_mRNA
        a20_m,                                     # A20_protein
        f_a20_c,                                   # A20_complex
        1 if (stim or tnf_s) else 0,               # TNFR_drive
        tnfr_d,                                    # TNFR1_active
        tnfr1,                                     # TRADD_recruited
        1 if (tradd and not a20_c) else 0,         # RIPK1_signal
        ripk1,                                     # TAK1_active
        tak1,                                      # IKK_active
        f_ikb_rec,                                 # IkB_recovery
        f_ikb,                                     # IkB_available
        1 if (ikk and not ikb) else 0,             # NFkB_released
        1 if (nf_rel and imp and not ikb) else 0,  # NFkB_nuclear
        1 if (nf_nuc and not ikb) else 0,          # NFkB_active
        nf_act,                                    # NFKBIA_mRNA
        nf_act,                                    # TNF_mRNA
        tnf_m,                                     # TNF_secreted
    )


def simulate_gillespie_custom(
    init_state_dict,
    max_time: float = 100.0,
    runs: int = 5000,
    seed: int = 1,
    a20_rule: str = 'strict_and',
    ikb_mode: str = 'D',
    rate_up=None,
    rate_down=None,
    itch_on_bias: float = 0.0,
    itch_off_bias: float = 0.0,
    a20_basal_leak: float = 0.0,
    ikb_resynth_leak: float = 0.0,
):
    if rate_up is None:
        r_up = [1.0] * 23
    elif isinstance(rate_up, (int, float)):
        r_up = [float(rate_up)] * 23
    else:
        r_up = [float(rate_up.get(n, 1.0)) for n in NODES]

    if rate_down is None:
        r_dn = [1.0] * 23
    elif isinstance(rate_down, (int, float)):
        r_dn = [float(rate_down)] * 23
    else:
        r_dn = [float(rate_down.get(n, 1.0)) for n in NODES]

    init_list = [int(init_state_dict[n]) for n in NODES]
    rng = random.Random(seed)
    res_count = 0
    nfkb_count = 0

    i_itch, i_a20p, i_a20c = 5, 7, 8
    i_tax, i_rnf, i_ikk, i_ikb = 2, 3, 14, 16
    i_nfa, i_ikbm, i_tnfm, i_tnfs = 19, 20, 21, 22

    for _ in range(runs):
        s = init_list[:]
        t = 0.0
        while t < max_time:
            f = eval_targets(s, a20_rule, ikb_mode)
            total = 0.0
            rates = [0.0] * 23
            for i in range(5, 23):
                si = s[i]
                fi = f[i]
                q = 0.0
                if si == 0 and fi == 1:
                    q += r_up[i]
                elif si == 1 and fi == 0:
                    q += r_dn[i]
                if i == i_itch:
                    if si == 0 and itch_on_bias > 0.0:
                        q += itch_on_bias
                    elif si == 1 and itch_off_bias > 0.0:
                        q += itch_off_bias
                elif i == i_a20c and a20_basal_leak > 0.0:
                    if si == 0 and s[i_a20p] and s[i_tax] and s[i_rnf]:
                        q += a20_basal_leak
                elif i == i_ikb and ikb_resynth_leak > 0.0:
                    if si == 0 and s[i_ikbm] and s[i_ikk]:
                        q += ikb_resynth_leak
                rates[i] = q
                total += q
            if total <= 0.0:
                break
            dt = -math.log(1.0 - rng.random()) / total
            if t + dt > max_time:
                break
            t += dt
            u = rng.random() * total
            acc = 0.0
            for i in range(5, 23):
                acc += rates[i]
                if acc >= u:
                    s[i] = 1 - s[i]
                    break
        if s[i_nfa] == 0 and s[i_ikk] == 0 and s[i_tnfm] == 0 and s[i_tnfs] == 0:
            res_count += 1
        if s[i_nfa] == 1:
            nfkb_count += 1
    return {
        'runs': runs,
        'resolved': res_count / runs,
        'resolved_ci95': wilson(res_count, runs),
        'nfkb_active': nfkb_count / runs,
        'nfkb_ci95': wilson(nfkb_count, runs),
    }


def established_init(itch_supply=0, rnf11=1, tax1bp1=1, nfkb_import=1):
    s = {n: 1 for n in NODES}
    s['STIMULUS'] = 0
    s['ITCH_supply'] = itch_supply
    s['TAX1BP1_available'] = tax1bp1
    s['RNF11_available'] = rnf11
    s['NFkB_import'] = nfkb_import
    s['ITCH_functional'] = 0
    s['IkB_available'] = 0
    s['IkB_recovery'] = 0
    s['A20_complex'] = 0
    return s


def main(runs: int = 5000):
    t0 = time.time()
    doses = [0.0, 0.01, 0.1, 1.0]
    results = {}
    seed_i = 1000

    # Part A: A20 complex rule variants
    part_a = []
    for label, rule, rnf, eps in [
        ('baseline_strict_AND_RNF11=1', 'strict_and', 1, 0.0),
        ('baseline_strict_AND_RNF11=0', 'strict_and', 0, 0.0),
        ('no_RNF11_required_RNF11=0',   'no_rnf11',   0, 0.0),
        ('scaffold_OR_RNF11=0_TAX=1',   'scaffold_or',0, 0.0),
        ('itch_or_scaffolds_RNF11=1',   'itch_or_scaffolds', 1, 0.0),
        ('partial_basal_A20_eps=0.005', 'strict_and', 1, 0.005),
        ('partial_basal_A20_eps=0.02',  'strict_and', 1, 0.02),
        ('partial_basal_A20_eps=0.05',  'strict_and', 1, 0.05),
    ]:
        for c in doses:
            seed_i += 1
            r = c / (1.0 + c)
            res = simulate_gillespie_custom(
                established_init(rnf11=rnf), max_time=100.0, runs=runs, seed=seed_i,
                a20_rule=rule, ikb_mode='D', itch_on_bias=r, a20_basal_leak=eps,
            )
            part_a.append({'variant': label, 'a20_rule': rule, 'rnf11': rnf, 'basal_leak': eps, 'c': c, 'r': r, **res})
    results['part_A_a20_rules'] = part_a

    # Part B: IkBa recovery modes
    part_b = []
    for label, mode, gam in [
        ('variant_D_baseline',        'D', 0.0),
        ('variant_D_resynth_g=0.01',  'D', 0.01),
        ('variant_D_resynth_g=0.05',  'D', 0.05),
        ('variant_D_resynth_g=0.10',  'D', 0.10),
        ('variant_M_induced_escapes', 'M_induced_escapes', 0.0),
        ('variant_R_full_recovery',   'R', 0.0),
    ]:
        for c in doses:
            seed_i += 1
            r = c / (1.0 + c)
            res = simulate_gillespie_custom(
                established_init(), max_time=100.0, runs=runs, seed=seed_i,
                a20_rule='strict_and', ikb_mode=mode, itch_on_bias=r, ikb_resynth_leak=gam,
            )
            part_b.append({'variant': label, 'ikb_mode': mode, 'resynth_leak': gam, 'c': c, 'r': r, **res})
    results['part_B_ikb_modes'] = part_b

    # Part C: Relative transition rates
    part_c = []
    fast_sig = {'TNFR_drive', 'TNFR1_active', 'TRADD_recruited', 'RIPK1_signal', 'TAK1_active',
                'IKK_active', 'A20_complex', 'NFkB_released', 'NFkB_nuclear', 'NFkB_active'}
    slow_tx = {'A20_mRNA', 'A20_protein', 'NFKBIA_mRNA', 'TNF_mRNA', 'TNF_secreted',
               'IkB_recovery', 'IkB_available'}
    fast5 = {n: (5.0 if n in fast_sig else 1.0) for n in NODES}
    slow02 = {n: (0.2 if n in slow_tx else 1.0) for n in NODES}
    for label, rup, rdn in [
        ('uniform_1_1',            1.0,    1.0),
        ('fast_signaling_5x',      fast5,  fast5),
        ('slow_transcription_0.2x',slow02, slow02),
        ('asymmetric_fast_off_2x', 1.0,    2.0),
        ('asymmetric_slow_off_0.5x',1.0,   0.5),
    ]:
        for c in doses:
            seed_i += 1
            r = c / (1.0 + c)
            res = simulate_gillespie_custom(
                established_init(), max_time=100.0, runs=runs, seed=seed_i,
                a20_rule='strict_and', ikb_mode='D', rate_up=rup, rate_down=rdn, itch_on_bias=r,
            )
            part_c.append({'variant': label, 'c': c, 'r': r, **res})
    results['part_C_rates'] = part_c

    # Part D: Equal stationary pi, different ITCH switching frequency
    part_d = []
    for pi_target, r0 in [(1 / 11, 0.1 / 1.1), (1 / 3, 0.5)]:
        for lam in [0.05, 0.2, 1.0, 5.0, 20.0]:
            seed_i += 1
            rdn_map = {n: 1.0 for n in NODES}
            rdn_map['ITCH_functional'] = lam
            r_on = lam * r0
            res = simulate_gillespie_custom(
                established_init(), max_time=100.0, runs=runs, seed=seed_i,
                a20_rule='strict_and', ikb_mode='D', rate_up=1.0, rate_down=rdn_map, itch_on_bias=r_on,
            )
            part_d.append({'pi': pi_target, 'lambda': lam, 'r_on': r_on, 'r_off': lam, **res})
    results['part_D_switching_freq'] = part_d

    out_path = Path(__file__).resolve().parent / "sensitivity_results.json"
    out_path.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f"Completed 86 conditions x {runs} trajectories in {time.time() - t0:.1f}s -> {out_path}")


if __name__ == '__main__':
    main()
