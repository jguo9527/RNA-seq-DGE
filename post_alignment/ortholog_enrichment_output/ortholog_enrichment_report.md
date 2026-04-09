# Ortholog-Based Enrichment for limma Interaction Genes

This report maps hemp interaction genes to Arabidopsis ortholog symbols using local annotation text plus MyGene lookup, then runs g:Profiler enrichment in Arabidopsis space.
The enrichment uses Arabidopsis genome-wide background, not a custom hemp transcriptome universe, so the results should be treated as pathway-level interpretation rather than exact odds-ratio calibration.

## Mapping summary

- Significant limma interaction genes: 698
- Interaction genes mapped to Arabidopsis-like orthologs: 209
- Unique Arabidopsis ortholog symbols from all interaction genes: 192
- Flowering-focused interaction genes mapped: 22
- Unique Arabidopsis ortholog symbols from flowering-focused interaction genes: 21

## Top enriched terms from all interaction genes

| source   | native     | name                                            |     p_value |   intersection_size |   term_size |
|:---------|:-----------|:------------------------------------------------|------------:|--------------------:|------------:|
| GO:BP    | GO:0022900 | electron transport chain                        | 6.14409e-10 |                  15 |         193 |
| GO:BP    | GO:0006091 | generation of precursor metabolites and energy  | 2.29662e-06 |                  17 |         462 |
| GO:BP    | GO:0006796 | phosphate-containing compound metabolic process | 0.00256873  |                  30 |        2025 |
| GO:BP    | GO:0006793 | phosphorus metabolic process                    | 0.00264646  |                  30 |        2028 |
| GO:BP    | GO:0006810 | transport                                       | 0.00268797  |                  34 |        2477 |
| GO:BP    | GO:0009987 | cellular process                                | 0.00271153  |                 122 |       16591 |
| GO:BP    | GO:0006119 | oxidative phosphorylation                       | 0.00288427  |                   5 |          40 |
| GO:BP    | GO:1902600 | proton transmembrane transport                  | 0.00336108  |                   8 |         155 |
| GO:BP    | GO:0051234 | establishment of localization                   | 0.00450165  |                  34 |        2536 |
| WP       | WP:WP2108  | Flower development initiation                   | 0.00494621  |                   3 |          17 |
| GO:BP    | GO:0022904 | respiratory electron transport chain            | 0.00721729  |                   5 |          48 |
| GO:BP    | GO:0051179 | localization                                    | 0.0103274   |                  34 |        2636 |

## Top enriched terms from flowering-focused interaction genes

| source   | native     | name                                                |     p_value |   intersection_size |   term_size |
|:---------|:-----------|:----------------------------------------------------|------------:|--------------------:|------------:|
| GO:BP    | GO:0032774 | RNA biosynthetic process                            | 3.17908e-05 |                  14 |        3395 |
| GO:BP    | GO:0141187 | nucleic acid biosynthetic process                   | 4.22386e-05 |                  14 |        3470 |
| GO:BP    | GO:0034654 | nucleobase-containing compound biosynthetic process | 9.36464e-05 |                  14 |        3690 |
| GO:BP    | GO:0016070 | RNA metabolic process                               | 0.000111792 |                  14 |        3741 |
| GO:BP    | GO:0006139 | nucleobase-containing compound metabolic process    | 0.000173268 |                  15 |        4625 |
| GO:BP    | GO:0008380 | RNA splicing                                        | 0.00025831  |                   6 |         355 |
| GO:BP    | GO:0090304 | nucleic acid metabolic process                      | 0.00055212  |                  14 |        4238 |
| GO:BP    | GO:0009639 | response to red or far red light                    | 0.000994672 |                   5 |         242 |
| GO:BP    | GO:0010467 | gene expression                                     | 0.00107273  |                  15 |        5284 |
| KEGG     | KEGG:03040 | Spliceosome                                         | 0.00116247  |                   5 |         205 |
| GO:BP    | GO:0006397 | mRNA processing                                     | 0.00128686  |                   6 |         468 |
| GO:BP    | GO:0009416 | response to light stimulus                          | 0.00147829  |                   7 |         764 |

## Interpretation

1. If plastid, thylakoid, and photosynthesis terms dominate the full interaction set, that indicates flowering stage differences are strongly coupled to cultivar-specific chloroplast and energy-state remodeling.
2. If circadian rhythm, developmental timing, hormone signaling, or RNA-splicing terms appear in the flowering-focused subset, those are the strongest candidates for causal flowering regulation rather than downstream metabolic state changes.
3. Use the flowering-focused mapped gene table first for breeding prioritization, and treat the broader chloroplast-enriched set as context for physiology and source-sink transition during floral induction.