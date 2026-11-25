# S1C–ITCH Inflammatory Response Modeling — Supplementary Data and Code

Supplementary data, structural coordinates, reconstructed gene network tables, 23-variable Boolean model specification, and sensitivity analysis code accompanying the study:

> **V.A. Ivanisenko, A.A. Nesterov, A.S. Venzel, P.S. Demenkov, T.V. Ivanisenko, V.N. Danilenko, N.A. Kolchanov.**  
> *Computational modeling of a potential mechanism underlying the anti-inflammatory action of Limosilactobacillus fermentum U-21 S1C protease through protein–protein interaction with the ubiquitin ligase ITCH.*  
> Bulletin of Biotechnology and Physicochemical Biology named after Yu.A. Ovchinnikov (*Вестник биотехнологии и физико-химической биологии им. Ю.А. Овчинникова*), 2026; 22(Special Issue).

---

## Repository Contents / Состав репозитория

| File | Description (EN) | Описание (RU) |
|---|---|---|
| `S1C_ITCH_complex_AF_Multimer_v3_rank001.pdb` | Atomic coordinates of the predicted S1C (`WP_003682316.1`, residues 1–423) – human ITCH N-lobe HECT fragment (`Q96J02`, residues 548–777) complex (`AlphaFold-Multimer v3`, unrelaxed `rank_001 model_5 seed_000`) | Атомные координаты предсказанного комплекса S1C–ITCH (`AlphaFold-Multimer v3`, `rank_001`) |
| `Table_S1_Candidate_pairs_ipSAE02.tsv` | Complete structural screening table of all 17 candidate pairs with `ipSAE ≥ 0.2` (16 with negative PyRosetta `dG_separated` + 1 without Rosetta score), threshold membership (`≥0.2`, `≥0.3`, `≥0.4`, `≥0.5`), PRIDE proteomics support, STRING interactors, and GO network membership | Полная таблица всех 17 кандидатных пар (`ipSAE ≥ 0,2`), оценки `dG_separated`, пороги `0,2/0,3/0,4/0,5`, поддержка PRIDE и вхождение в сеть GO |
| `Table_S2_S1C_ITCH_interface_contacts_4A.tsv` | All 33 interchain residue pairs within 4.0 Å heavy-atom distance in the predicted S1C–ITCH complex, highlighting the 6 N-lobe residues shared with the autoinhibitory 2,3-linker interface in PDB `5XMC` (`Arg594`, `Glu595`, `Phe598`, `His602`, `Phe654`, `Asp656`) | Перечень всех 33 контактных пар остатков S1C–ITCH (≤ 4,0 Å) и 6 общих остатков с аутоингибиторным интерфейсом PDB `5XMC` |
| `Table_S3_Network_nodes_685.tsv` | Complete list of all 685 nodes in the reconstructed ANDSystem inflammatory response network (`dbversion 20260417`) with degree, normalized centrality `C_D`, permutation `p`-value (5,000 `G(n,m)` graphs, `seed = 517`), `q_BH`, hub status, 27-node core membership, and Boolean model inclusion | Полный перечень 685 вершин реконструированной сети ANDSystem (`dbversion 20260417`) со степенями, центральностью, `p`, `q_BH` и вхождением в ядро |
| `Table_S4_Network_edges_7387.tsv` | Complete list of all 7,387 directed typed regulatory edges with original ANDSystem relation IDs and signs | Полный перечень 7387 направленных типизированных связей сети с исходными идентификаторами ANDSystem |
| `Table_S5_Regulatory_core_27_nodes_283_edges.tsv` | Complete list of all 283 directed edges in the 27-node regulatory core (`NetworkX 3.6.1`, Louvain `resolution = 1`, `weight = None`, consensus of `seed = 1, 11, 23` + restored `TNF`) | Перечень 283 направленных связей 27-вершинного регуляторного ядра |
| `Table_S6_Boolean_model_23_variables_rules_initial_states.tsv` | Full specification of all 23 Boolean model variables, Variant D and Variant R update functions, initial states (established-active and resting), and primary literature provenance | Полная спецификация 23 переменных, логических функций (варианты D и R) и начальных состояний булевой модели |
| `Table_S7_Boolean_sensitivity_analysis.tsv` | Full results of the 86-condition × 5,000-trajectory sensitivity analysis (A20 complex rule variants, IκBα recovery modes, relative transition rates / timescale separation, and ITCH switching frequency at fixed stationary fraction `π`) with 95% Wilson CIs | Полные результаты анализа чувствительности (86 условий × 5000 траекторий) с 95%-ными ДИ Уилсона |
| `sensitivity_analysis.py` | Standalone Python 3 (standard library only) Gillespie simulation script reproducing all 86 sensitivity conditions | Автономный скрипт на Python 3 (только стандартная библиотека) для воспроизведения всех 86 условий анализа чувствительности |
| `sensitivity_results.json` | Machine-readable JSON output of `sensitivity_analysis.py` | Машиночитаемые результаты расчёта анализа чувствительности |
| `README_Воспроизводимость.txt` | Russian reproducibility manifest with software/database versions and SHA-256 checksums | Манифест воспроизводимости на русском языке с контрольными суммами SHA-256 |

---

## Reproducing the Sensitivity Analysis / Воспроизведение расчётов

Requires Python 3.9+ (no third-party packages required):

```bash
python3 sensitivity_analysis.py
```

The script runs all 86 conditions (5,000 trajectories each, deterministic seeds) and writes `sensitivity_results.json` (~2 minutes on a modern laptop).

---

## External Data Sources / Внешние источники данных

- **PRIDE Proteomics:** [PXD050857](https://www.ebi.ac.uk/pride/archive/projects/PXD050857) (*L. fermentum* U-21 extracellular vesicles, `proteinGroups453.xlsx`)
- **NCBI RefSeq Genome:** [NZ_CP103293.1](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP103293.1) (*L. fermentum* U-21; S1C protease `WP_003682316.1`, locus `C0965_RS00130`)
- **UniProtKB:** [Q96J02](https://www.uniprot.org/uniprotkb/Q96J02/entry) (human E3 ubiquitin-protein ligase ITCH)
- **RCSB PDB:** [5XMC](https://www.rcsb.org/structure/5XMC) (autoinhibited ITCH), [3JW0](https://www.rcsb.org/structure/3JW0) (NEDD4L HECT–UbcH5B complex)

---

## License

- Code (`sensitivity_analysis.py`): **MIT License** (see `LICENSE`).
- Supplementary data tables and coordinate files: **CC BY 4.0**.
