# Hemp Flowering Enrichment Summary With Term Explanations

## What this report covers

This report explains the enriched pathway terms from the limma-confirmed interaction genes after ortholog mapping to Arabidopsis symbols.
It summarizes two layers:

1. All significant interaction genes (broad cultivar-specific stage response)
2. Flowering-focused interaction genes (candidate flowering regulators)

Use the flowering-focused section for breeding target prioritization, and the all-interactions section for whole-plant physiological context.

## Dataset context

- Significant limma interaction genes: 698
- Mapped interaction genes for enrichment: 209 (192 unique Arabidopsis-like symbols)
- Flowering-focused mapped genes: 22 (21 unique symbols)

---

## A. All Interaction Genes: Term-by-Term Explanations

| Term | ID | Plain-language explanation | Why it matters for hemp flowering biology |
|---|---|---|---|
| electron transport chain | GO:0022900 | A sequence of protein complexes that move electrons and generate cellular energy through redox reactions. | Strong enrichment suggests cultivar-dependent differences in energy conversion capacity during flowering transition. |
| generation of precursor metabolites and energy | GO:0006091 | Core metabolism that creates ATP and molecular building blocks needed for growth and development. | Indicates stage responses include shifts in energy allocation and metabolic readiness for reproductive development. |
| phosphate-containing compound metabolic process | GO:0006796 | Reactions involving phosphate-containing molecules such as ATP, nucleotides, and phosphorylated intermediates. | Highlights broad signaling and energy-transfer rewiring among cultivars. |
| phosphorus metabolic process | GO:0006793 | A broad category for metabolic reactions involving phosphorus compounds. | Supports large-scale shifts in phosphorylation and ATP-dependent regulation. |
| transport | GO:0006810 | Movement of ions, metabolites, and macromolecules across membranes or between compartments. | Suggests flowering stage differences involve resource redistribution and signaling transport. |
| cellular process | GO:0009987 | A high-level umbrella term for cell-level biological activities. | Confirms interaction effects are widespread, not limited to a single pathway branch. |
| oxidative phosphorylation | GO:0006119 | ATP production linked to electron transport and proton gradients. | Points to cultivar-specific differences in mitochondrial/chloroplast energy efficiency during floral transition. |
| proton transmembrane transport | GO:1902600 | Movement of protons across membranes, which powers ATP synthesis. | Mechanistically tied to bioenergetics and chloroplast/mitochondrial state changes. |
| establishment of localization | GO:0051234 | Processes that place molecules or organelles at specific cellular locations. | Implies dynamic subcellular reorganization during stage transition differs by cultivar. |
| Flower development initiation | WP:WP2108 | Pathway-level gene set associated with the onset of floral development. | Direct evidence that interaction genes include true flowering-initiation components. |
| respiratory electron transport chain | GO:0022904 | Electron transfer steps specifically linked to respiratory metabolism. | Supports cultivar-dependent respiratory tuning as plants move into reproductive stage. |
| localization | GO:0051179 | General placement or maintenance of molecules and structures in specific locations. | Adds evidence for cellular architecture and trafficking changes around flowering. |

### Interpretation of section A

The broad interaction set is dominated by energy and transport programs, meaning flowering stage effects are not only developmental-gene effects but also whole-system physiological remodeling that depends on genetic background.

---

## B. Flowering-Focused Interaction Genes: Term-by-Term Explanations

| Term | ID | Plain-language explanation | Why it matters for hemp flowering biology |
|---|---|---|---|
| RNA biosynthetic process | GO:0032774 | Production of RNA molecules, including transcription and related synthesis routes. | Indicates cultivar-specific flowering responses strongly involve transcriptional control points. |
| nucleic acid biosynthetic process | GO:0141187 | Synthesis of DNA/RNA-related nucleic acid molecules. | Suggests stage transition engages regulatory programs at nucleic acid production level. |
| nucleobase-containing compound biosynthetic process | GO:0034654 | Formation of nucleotides, nucleic acids, and related compounds. | Reinforces active nucleic-acid synthesis and gene-expression infrastructure shifts. |
| RNA metabolic process | GO:0016070 | Processing, modification, turnover, and handling of RNA molecules. | Points to post-transcriptional regulation as a key driver of cultivar-specific flowering behavior. |
| nucleobase-containing compound metabolic process | GO:0006139 | Broad metabolism of nucleic-acid-related compounds. | Consistent with high regulatory demand during developmental switching. |
| RNA splicing | GO:0008380 | Removal of introns and joining of exons to produce mature RNAs. | Mechanistic signal that isoform-level regulation may shape flowering timing and plasticity. |
| nucleic acid metabolic process | GO:0090304 | All metabolism related to DNA/RNA molecules. | Confirms nucleic-acid regulatory machinery is central in interaction effects. |
| response to red or far red light | GO:0009639 | Cellular or organismal response to red/far-red wavelengths, often phytochrome-mediated. | Directly relevant to photoperiod sensing and flowering induction logic. |
| gene expression | GO:0010467 | Overall process from transcription to mature RNA/protein output. | Indicates that transcriptional program timing differs across cultivar-stage combinations. |
| Spliceosome | KEGG:03040 | Core RNA-protein machinery that carries out splicing reactions. | Strong pathway evidence for alternative splicing as a flowering-regulation lever. |
| mRNA processing | GO:0006397 | RNA maturation steps including capping, splicing, and transcript processing. | Supports regulation occurring between transcription and protein synthesis. |
| response to light stimulus | GO:0009416 | Any biological response triggered by light. | Connects interaction signal to photoreceptor and light-driven flowering pathways. |

### Interpretation of section B

The flowering-focused subset is enriched for two tightly linked regulatory axes:

1. Light and photoperiod signal perception/response
2. RNA processing and splicing control of gene output

This supports a model in which cultivar-specific flowering behavior is controlled by how light cues are interpreted and translated into transcript-isoform and gene-expression programs.

---

## Practical implications for variety improvement

1. Prioritize markers in light-response and RNA-splicing regulators because these terms are enriched in the flowering-focused subset rather than only in broad metabolism.
2. Treat energy/transport-enriched genes as supporting physiological context markers; they may correlate with flowering state but are less likely to be primary flowering switches.
3. Build selection panels around ortholog-backed candidates appearing in these enriched terms, especially genes aligned with CAULIFLOWER, AGL16, SVP, ELF3, CRY1, SPL12, GA20OX1, CYP707A2, ARF/IAA/SAUR, and SR/RS splicing factors.

## Caveat

Term enrichment was performed in Arabidopsis ortholog space and uses Arabidopsis background. This is ideal for mechanistic interpretation but should be combined with hemp-specific validation before final deployment in breeding decisions.