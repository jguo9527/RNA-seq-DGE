# Hemp Flowering limma-voom Confirmation

This run confirms the hemp flowering marker set with edgeR filtering, TMM normalization, and limma-voom with quality weights using the same stage plus cultivar interaction design as the Python analysis.

## Dataset summary

- Samples analyzed: 30
- Genes retained after filterByExpr: 18802
- Shared stage genes at FDR < 0.05: 9649
- Interaction genes at FDR < 0.05: 698
- SGxAV correlation to SG: 0.588
- SGxAV correlation to EAV: 0.405

Mean QC by cultivar and stage:

| cultivar | stage | genic_fraction | mapped_fraction |
|---|---|---|---|
| EAV | before | 0.556 | 0.864 |
| EAV | after | 0.452 | 0.832 |
| LAV | before | 0.478 | 0.832 |
| LAV | after | 0.534 | 0.838 |
| Lifter | before | 0.363 | 0.785 |
| Lifter | after | 0.860 | 0.952 |
| SG | before | 0.530 | 0.831 |
| SG | after | 0.507 | 0.863 |
| SGxAV | before | 0.635 | 0.875 |
| SGxAV | after | 0.456 | 0.825 |

## Cross-method confirmation

- Python vs limma shared-stage overlap at FDR < 0.05: 7025 genes
- Python vs limma interaction overlap at FDR < 0.05: 248 genes
- Overlap of top 100 shared-stage genes: 70
- Overlap of top 100 interaction genes: 40

## Top shared flowering markers

| gene_id | gene | description | logFC | adj.P.Val |
|---|---|---|---|---|
| LOC115705224 | LOC115705224 | protein LATE ELONGATED HYPOCOTYL | -4.506490 | 1.186308e-22 |
| LOC115714162 | LOC115714162 | adagio protein 3 |  3.291159 | 5.449893e-22 |
| LOC115711591 | LOC115711591 | two-component response regulator-like APRR1 |  1.876239 | 2.006795e-19 |
| LOC115714579 | LOC115714579 | B-box zinc finger protein 18 | -3.177430 | 1.658521e-18 |
| LOC115702321 | LOC115702321 | cold-regulated protein 27 |  4.507792 | 3.023338e-18 |
| LOC133029979 | LOC133029979 | zinc finger protein CONSTANS-LIKE 15-like |  2.914468 | 3.226140e-18 |
| LOC115721095 | LOC115721095 | lysine-specific demethylase JMJ30 |  2.104978 | 1.115508e-17 |
| LOC115714985 | LOC115714985 | ABC transporter G family member 5 | -1.901940 | 4.012369e-17 |
| LOC115712781 | LOC115712781 | REF/SRPP-like protein At1g67360 | -1.809875 | 7.222137e-17 |
| LOC115706039 | LOC115706039 | protein EARLY FLOWERING 4 |  3.825627 | 1.140852e-16 |
| LOC115705704 | LOC115705704 | protein LNK3 | -4.000734 | 2.773077e-16 |
| LOC115718944 | LOC115718944 | transcription factor BOA |  2.471268 | 3.026122e-16 |
| LOC115694953 | LOC115694953 | probable anion transporter 6, chloroplastic | -1.846255 | 3.769073e-16 |
| LOC115699196 | LOC115699196 | alpha-amylase 3, chloroplastic |  1.684275 | 4.343457e-16 |
| LOC115721126 | LOC115721126 | uncharacterized LOC115721126 |  3.542277 | 6.026888e-16 |

## Top cultivar-specific interaction genes

| gene_id | gene | description | F | adj.P.Val | effect_span |
|---|---|---|---|---|---|
| TRNAS-GGA | TRNAS-GGA | transfer RNA serine (anticodon GGA) | 29.09435 | 0.0001278998 | 0.4802993 |
| LOC133030118 | LOC133030118 | acetyl-coenzyme A carboxylase carboxyl transferase subunit beta, chloroplastic-like | 26.15655 | 0.0001743406 | 1.2862199 |
| TRNAG-GCC_10 | TRNAG-GCC | transfer RNA glycine (anticodon GCC) | 25.50262 | 0.0001743406 | 0.7012720 |
| LOC133037045 | LOC133037045 | uncharacterized LOC133037045 | 24.74677 | 0.0001743406 | 0.5581238 |
| TRNAF-GAA_3 | TRNAF-GAA | transfer RNA phenylalanine (anticodon GAA) | 22.98458 | 0.0001743406 | 0.5276304 |
| LOC133035721 | LOC133035721 | cytochrome f-like | 23.74304 | 0.0001743406 | 0.3710928 |
| TRNAC-GCA_6 | TRNAC-GCA | transfer RNA cysteine (anticodon GCA) | 24.54046 | 0.0001743406 | 0.3466375 |
| LOC133037471 | LOC133037471 | photosystem I P700 chlorophyll a apoprotein A1 | 24.36034 | 0.0001777040 | 0.3896795 |
| LOC115700311 | LOC115700311 | DNA-directed RNA polymerase subunit alpha-like | 23.45617 | 0.0001919564 | 0.2939548 |
| LOC133037037 | LOC133037037 | photosystem II CP47 reaction center protein-like | 23.20054 | 0.0002369811 | 0.2725917 |
| LOC133037460 | LOC133037460 | NAD(P)H-quinone oxidoreductase subunit 6, chloroplastic | 22.37107 | 0.0002563705 | 0.4767492 |
| LOC133037468 | LOC133037468 | cytochrome b6 | 20.11736 | 0.0003681438 | 1.7310754 |
| A5N79_gp09 | rpl2 | NA | 21.27670 | 0.0003681438 | 0.3760484 |
| LOC133035722 | LOC133035722 | chloroplast envelope membrane protein | 21.12716 | 0.0003720878 | 0.5822867 |
| LOC133037474 | LOC133037474 | ATP synthase subunit beta, chloroplastic | 20.78013 | 0.0003796603 | 0.7075683 |
| LOC133037044 | LOC133037044 | uncharacterized LOC133037044 | 20.12644 | 0.0004481961 | 0.4770749 |
| LOC133039687 | LOC133039687 | photosystem I P700 chlorophyll a apoprotein A2-like | 18.68828 | 0.0004768807 | 2.6525426 |
| A5N79_gp35 | atp1 | NA | 19.87875 | 0.0004954689 | 0.3504865 |
| LOC133037473 | LOC133037473 | ATP synthase subunit alpha, chloroplastic | 18.60850 | 0.0005037553 | 2.1123284 |
| LOC115709109 | LOC115709109 | heavy metal-associated isoprenylated plant protein 7 | 18.25388 | 0.0005037553 | 1.4654938 |

## Flowering-focused candidate genes

| gene_id | gene | description | EAV_log2_fc | SG_log2_fc | SGxAV_log2_fc | interaction_fdr |
|---|---|---|---|---|---|---|
| LOC115702300 | LOC115702300 | truncated transcription factor CAULIFLOWER A |  0.82607488 |  9.65661185 |  9.89163914 | 0.001334810 |
| LOC115702792 | LOC115702792 | serine/arginine-rich splicing factor SR45a |  0.07719351 | -0.78003137 | -0.57845121 | 0.004264373 |
| LOC115695700 | LOC115695700 | protein transport protein SFT2 |  0.21674337 | -0.29880322 | -0.05834071 | 0.010084963 |
| LOC115725225 | LOC115725225 | auxin-responsive protein SAUR21 |  0.38247236 |  0.02403190 |  0.09954373 | 0.010084963 |
| LOC115718499 | LOC115718499 | abscisic acid 8'-hydroxylase 2 | -0.31566972 | -1.13491649 |  0.02787571 | 0.010482875 |
| LOC115711648 | LOC115711648 | zinc finger protein CONSTANS-LIKE 6 |  0.29538942 | -0.36847611 |  0.71216536 | 0.020394445 |
| LOC115719221 | LOC115719221 | pre-mRNA-splicing factor 38 | -0.49029334 | -0.03865757 | -0.15111352 | 0.023609324 |
| LOC115697537 | LOC115697537 | alpha,alpha-trehalose-phosphate synthase [UDP-forming] 1 |  0.20031591 | -1.06234083 | -0.52156269 | 0.023912242 |
| LOC115704768 | LOC115704768 | auxin-responsive protein IAA27 | -0.47768867 | -0.17153970 |  0.34931300 | 0.029075947 |
| LOC115718181 | LOC115718181 | serine/arginine-rich splicing factor RS40 | -0.89746066 | -0.18947166 | -0.37972958 | 0.035080132 |
| LOC115721274 | LOC115721274 | MADS-box protein SVP | -0.19472839 | -0.70396578 | -0.92705624 | 0.042606340 |
| LOC115722191 | LOC115722191 | pre-mRNA splicing factor SR-like 1 | -0.29590394 | -0.11315376 | -0.11966688 | 0.060104960 |
| LOC115718979 | LOC115718979 | auxin-responsive protein IAA16 | -0.07848380 |  0.01719951 | -0.03944414 | 0.069006866 |
| LOC115712114 | LOC115712114 | gibberellin 20 oxidase 1-B | -0.76042074 | -0.86787796 | -1.62004132 | 0.069769887 |
| LOC115713961 | LOC115713961 | serine/arginine-rich splicing factor RSZ21 |  0.62433703 | -0.21379306 |  0.12923620 | 0.070884134 |
| LOC115698037 | LOC115698037 | probable pre-mRNA-splicing factor ATP-dependent RNA helicase DEAH2 | -0.53285951 | -0.08073480 | -0.10362675 | 0.074955563 |
| LOC115703149 | LOC115703149 | protein EARLY FLOWERING 3 |  0.13551564 |  0.47859346 |  0.42811958 | 0.076315267 |
| LOC115705617 | LOC115705617 | auxin response factor 19 | -0.44531130 |  0.07863116 |  0.25766954 | 0.081504780 |
| LOC115709034 | LOC115709034 | agamous-like MADS-box protein AGL16 | -0.73813733 |  4.45179601 |  4.34072949 | 0.081770001 |
| LOC115705305 | LOC115705305 | serine/arginine-rich splicing factor RS2Z32 | -1.01887867 | -0.34372674 | -0.33536376 | 0.086485253 |
| LOC115705343 | LOC115705343 | auxin-responsive protein IAA17 |  0.36123394 |  0.16594115 | -0.40737871 | 0.087222446 |
| LOC115708788 | LOC115708788 | auxin-responsive protein SAUR71 |  0.30970135 |  0.89778917 |  0.23186829 | 0.087222446 |
| LOC115695743 | LOC115695743 | phytochrome A-associated F-box protein | -1.91105870 | -2.53111912 | -1.38610529 | 0.090273850 |
| LOC115717098 | LOC115717098 | probable trehalose-phosphate phosphatase F | -1.28904673 | -0.41623272 | -0.46665404 | 0.090427657 |
| LOC115709153 | LOC115709153 | abscisic acid 8'-hydroxylase CYP707A2 | -1.01169145 |  1.01951935 |  0.20403310 | 0.091850122 |

