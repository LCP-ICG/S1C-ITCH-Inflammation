ДОПОЛНИТЕЛЬНЫЕ МАТЕРИАЛЫ И АРХИВ ВОСПРОИЗВОДИМОСТИ
К статье В.А. Иванисенко и соавт. (Вестник биотехнологии и физико-химической биологии им. Ю.А. Овчинникова, 2026)
Дата формирования: 25.11.2026

СОСТАВ МАТЕРИАЛОВ:
1. S1C_ITCH_complex_AF_Multimer_v3_rank001.pdb — атомные координаты предсказанной модели комплекса S1C (WP_003682316.1, ост. 1–423) и фрагмента N-доли HECT-домена ITCH человека (UniProt Q96J02, ост. 548–777), AlphaFold-Multimer v3 (unrelaxed rank_001 model_5 seed_000).
2. Table_S1_Candidate_pairs_ipSAE02.tsv — полные результаты отбора всех 17 пар с ipSAE >= 0,2 (16 пар с отрицательной расчётной энергией dG_separated в PyRosetta ref2015 и 1 пара WP_023466532.1–CYP1A1 с ipSAE = 0,6685 без оценки Rosetta), с индикаторами прохождения порогов ipSAE >= 0,2; 0,3; 0,4; 0,5, протеомной поддержкой PRIDE PXD050857 и присутствием в очищенной сети воспалительного ответа.
3. Table_S2_S1C_ITCH_interface_contacts_4A.tsv — перечень всех 33 контактных пар остатков S1C–ITCH (расстояние между тяжёлыми атомами <= 4,0 Å), выделение 6 общих остатков с аутоингибиторным интерфейсом WW2/2,3-линкера ITCH (PDB 5XMC).
4. Table_S3_Network_nodes_685.tsv — полный перечень 685 вершин ориентированной сети регуляции воспалительного ответа ANDSystem (dbversion 20260417) с показателями степени, центральности C_D, p-значениями (5000 случайных графов G(n,m), seed=517), q_BH и принадлежностью к 27-вершинному регуляторному ядру.
5. Table_S4_Network_edges_7387.tsv — полный перечень 7387 направленных типизированных связей сети с исходными идентификаторами взаимодействий ANDSystem.
6. Table_S5_Regulatory_core_27_nodes_283_edges.tsv — перечень 283 направленных связей 27-вершинного регуляторного ядра (консенсус Лувена в NetworkX 3.6.1, resolution=1, weight=None, seed=1, 11, 23 + восстановление TNF).
7. Table_S6_Boolean_model_23_variables_rules_initial_states.tsv — полная спецификация 23 переменных, логических функций (варианты D и R) и начальных состояний булевой модели.
8. Table_S7_Boolean_sensitivity_analysis.tsv — полные результаты анализа чувствительности булевой модели (86 условий × 5000 траекторий: варианты правила комплекса A20, режимы восстановления IκBα, относительные интенсивности переходов и разделение временных масштабов, частота переключений ITCH при фиксированной стационарной доле π).
9. sensitivity_analysis.py и sensitivity_results.json — исходный скрипт и машиночитаемые результаты анализа чувствительности (Python 3.13, кросс-валидация с Pathway_boolean_modelling v0.14.1).

ВЕРСИИ БАЗ ДАННЫХ И ПРОГРАММНОГО ОБЕСПЕЧЕНИЯ:
- NCBI RefSeq / геном Limosilactobacillus fermentum U-21: NZ_CP103293.1 (проверка 11.10.2026)
- PRIDE proteomics archive: PXD050857 (proteinGroups453.xlsx, mqpar.xml)
- UniRef30_2020_06, HH-suite3 (HHblits), ColabFold (AlphaFold2-ptm, AlphaFold-Multimer v3)
- AlphaFold Protein Structure Database (модели доменов AF-...-model_v6; протокол отбора доменов по TED pLDDT >= 70)
- MaSIF-search / dMaSIF; ipSAE (Dunbrack, 2025); PyRosetta 4 (FastRelax, ref2015, InterfaceAnalyzer dG_separated)
- PyMOL 3.1.0; PDB 5XMC (аутоингибированный ITCH), PDB 3JW0 (NEDD4L–UbcH5B)
- QuickGO / Gene Ontology (GO:0006954, GO:0050727, выгрузка 10.10.2026, taxon 9606, Swiss-Prot, без IEA/NOT/ND/RCA)
- ANDSystem (формат v3, база знаний dbversion 20260417)
- Pathway_analysis, Pathway_boolean_modelling v0.14.1, NetworkX 3.6.1

КОНТРОЛЬНЫЕ СУММЫ SHA-256:
  263565299df769cc79dde81a8f5896aec607931fc82a12f4d00881611fb0dfea  LICENSE
  13f770f93c7265e1f66312166e7b36364178dc9454318a91bd80a9ef7de72c6f  README.md
  f9c8b8b675370270ed36966927c2fdaafdefde773daae9ec0f5e6c62617d91ee  S1C_ITCH_complex_AF_Multimer_v3_rank001.pdb
  92ccc3eed9ee131d8dec1573c7e0c55c3effbb3762726d9b66cec4ec99abd923  Table_S1_Candidate_pairs_ipSAE02.tsv
  6af53b0b3ab5b27429a62b2eed382ab4009163f5ee30e7d50ab7f88ee6e45f69  Table_S2_S1C_ITCH_interface_contacts_4A.tsv
  b621bda6d5eb4818528d5bef446a0c35c3dd0900cf9a942bcf63500806ff7ebb  Table_S3_Network_nodes_685.tsv
  52bda3b5974a3cd2658dbfcd736a1e7803b594c508193e5302cae250c9b1ca2d  Table_S4_Network_edges_7387.tsv
  c7ec4e4cbb50784c8adf28d5c613e608108fb70c275fd356c94da03e389cee16  Table_S5_Regulatory_core_27_nodes_283_edges.tsv
  1f9250bfcee00a81023cc73d5b549d2a1c7bd290611a99dccc20f975064e504d  Table_S6_Boolean_model_23_variables_rules_initial_states.tsv
  c92990ef33d4753a1de6692736979b35e00c635bd9b692788e5ce8fc70d00c40  Table_S7_Boolean_sensitivity_analysis.tsv
  60530867db7be822ff34b2a52c99847476db46624a8db3b00c01b69dc62c97cf  sensitivity_analysis.py
  7503046f4f46be77b703ebd6071ce9ce257bab8bb41402d41ed641826117d476  sensitivity_results.json
