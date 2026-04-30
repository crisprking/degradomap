# degradomap: methods note

**Author:** Abraham Trueba

**Date:** April 2026

## Abstract

We assessed whether features derivable from public databases (AlphaFold-predicted
structural druggability metrics, DepMap CRISPR gene-effect (Chronos) essentiality,
and cell-line expression breadth) can discriminate the 12 PROTAC-validated human
E3 ligases from 814 other human E3 candidate proteins. Using leave-one-out
cross-validated logistic regression, we observed AUCs of 0.473 (structural
features only), 0.587 (essentiality and expression), and 0.529 (combined). The
top-ranked candidates from the combined model were dominated by core ubiquitin
pathway components (SKP1, ELOB, DDB1, E2 conjugating enzymes) rather than novel
PROTAC E3 candidates. We conclude that the features captured by current public
structural and functional-genomics data are insufficient to predict PROTAC E3
tractability.

## Introduction

I came to this question through training in chemical biology (UC Berkeley) and bioinformatics. PROTAC E3 selection is the bottleneck of the field — the chemistry community has expanded the warhead toolkit dramatically, but we still rely on a handful of E3 ligases (CRBN, VHL) for almost all published degraders. I wanted to test whether public data — AlphaFold structures, DepMap functional genomics — could rank the unexploited E3s in a way that would help triage which ones are worth pursuing. The answer turned out to be more nuanced than yes or no, which is the subject of this and follow-up posts.

The human ubiquitin-proteasome system contains hundreds of E3 ligases, yet
PROTACs and molecular glues to date have engaged a small minority. Fewer than
20 proteins have published, validated degrader chemistry. The mechanistic
reasons are partly understood (existence of a tractable small-molecule binder,
substrate-adapter geometry, non-essential cellular role), but a quantitative,
data-driven approach to ranking the remaining ~600 human E3 candidates by
tractability has not been published. We tested whether public data alone can
provide such a ranking.

## Methods

### E3 ligase candidate set

We queried UniProt (REST API, accessed April 2026) for reviewed human entries
matching any of four query classes: (i) direct E3 ubiquitin ligases (GO:0061630,
GO:0004842, or "E3 ubiquitin"/"ubiquitin ligase" in function), (ii) substrate
adapters (DDB1/CUL4 association, DCAF protein names, GO:1990756), (iii) F-box
family proteins, (iv) BTB and SOCS box proteins. After deduplication, 859
candidates were obtained. Thirteen published PROTAC E3 ligases (CRBN, VHL, XIAP,
BIRC2, BIRC3, MDM2, DCAF15, DCAF16, DCAF1, RNF114, RNF4, KEAP1, FEM1B) were
resolved to UniProt accessions by gene-symbol lookup and added if missing.
Final candidate set: 872 proteins, 13 labeled positives.

### Structural features

AlphaFold-predicted structures were retrieved via the EMBL-EBI prediction API,
using the API to resolve current model versions rather than hardcoded URLs.
866/872 candidates had available predictions.

For each structure we computed per-residue solvent accessible surface area
(Shrake-Rupley algorithm via Biopython) and normalized to empirical maximum SASA
values from Tien et al. 2013. Analysis was restricted to residues with AlphaFold
pLDDT > 70. Ten features were computed: compactness (radius of gyration over
sqrt of confident-residue count), hydrophobic and aromatic surface fractions,
maximum hydrophobic surface patch size and count of patches of 5+ residues
(via connected components at 7 Angstrom threshold), pocket-like residue density
(partially-buried hydrophobic, 0.05 < relSASA < 0.40), maximum pocket cluster
size (8 Angstrom threshold), size-normalized version, plus mean pLDDT and
confident fraction.

### DepMap features

We used DepMap Public 26Q1 (released April 2026). For each candidate gene present
in the CRISPR gene-effect matrix (n=826 of 872 matched), we computed mean and
median Chronos score across all screened cell lines, plus fraction-essential at
thresholds Chronos < -0.5 and Chronos < -1.0. From the expression matrix
(n=860 matched) we computed mean log-TPM and fraction of cell lines with
expression above log2(TPM+1) > 1 and > 3.

### Discrimination experiment

826 proteins had complete data across all three feature classes (12 of 13
positives; DCAF16 was excluded because its high-confidence structural region
was below the 30-residue minimum). We trained logistic regression classifiers
(class-balanced, L2 regularization C=1) under leave-one-out cross-validation
across three feature sets: structural only (10 features), biological only
(4 features), combined (14 features).

## Results

### Per-feature comparison

Of 14 features, three showed statistically significant differences (p < 0.05)
between PROTAC-validated E3s and other candidates:

| Feature | Positives | Others | p |
|---|---|---|---|
| max_hp_patch | 6.50 | 10.52 | <0.001 |
| frac_expressed | 0.978 | 0.810 | <0.001 |
| mean_log_tpm | 4.67 | 3.55 | 0.001 |
| n_hp_patches_5plus | 1.92 | 2.83 | 0.015 |

Notably, the structural features that differed went in the OPPOSITE direction
from druggability theory: PROTAC-validated E3s have SMALLER and FEWER surface
hydrophobic patches than non-validated candidates. The biological features
behaved as predicted: PROTAC E3s are more broadly expressed across cell lines.

### Cross-validated discrimination

| Feature set | LOO-CV AUC |
|---|---|
| Structural (10 features) | 0.473 |
| Biological (4 features) | 0.587 |
| Combined (14 features) | 0.529 |

No model substantially exceeded random performance.

### Ranked predictions

The top-20 candidates ranked by the combined model consisted entirely of false
positives (no PROTAC-validated E3s in the top 80). The top-ranked proteins
were components of E3 ligase machinery rather than substrate adapters: SKP1,
ELOB, ELOC (cullin-RING adapters), DDB1 (CRBN/DCAF scaffold), UBE2D3, UBE2C,
UBE2V1, UBE2Z (E2 conjugating enzymes), ANAPC5, ANAPC11, CDC16 (APC/C complex
components).

### Bimodality of the PROTAC-validated set

The 12 positives split into two groups in the combined-model score
distribution: six (CRBN, MDM2, BIRC2, DCAF15, XIAP, VHL) ranked in the top
300, and six (DCAF1, RNF4, RNF114, KEAP1, BIRC3, FEM1B) ranked below 500.
This reflects the heterogeneity of the positive set: it contains both
classical pocket-druggable E3s (KEAP1, VHL, MDM2) and substrate-mimic-druggable
E3s (CRBN, DCAF15, where IMiDs and indisulam bind a tri-tryptophan cage rather
than a classical pocket). DCAF1 is an additional outlier with Chronos score
of -2.03 (highly essential), possibly reflecting that DCAF1's essential role
derives from its serine/threonine kinase activity rather than its
substrate-adapter function.

## Discussion

The most defensible interpretation of the v1 null result is that structural and functional-genomics features alone capture the protein side of PROTAC tractability but miss the chemistry side entirely. CRBN is not predictable from its structure or its DepMap profile; it is predictable from the historical fact that thalidomide was found to bind it. ChEMBL ligandability features (number of distinct scaffolds, best Kd, presence of co-crystal structures) are the obvious next addition.

A subsequent analysis (see follow-up posts) showed that splitting the 12 positives by binding mechanism — pocket binders, molecular glues, covalent — recovers signal that the pooled model erased. External validation against post-2024 PROTAC E3s confirmed the pocket-binder model partially generalizes, with a clear failure mode on Kelch β-propeller folds. These follow-up findings are documented in subsequent commits and writeups.

The principal limitation is the small positive set (n=12, mechanistically
heterogeneous). A second limitation is the absence of warhead-side features.
Tractability is a property of the protein-ligand pair, not of the protein
alone. Features like best-known small-molecule binder Kd or number of distinct
chemical scaffolds reported as binders might capture the missing signal but
require ChEMBL ligand-binding curation we did not perform here.

## References

Bekes, M., Langley, D. R., & Crews, C. M. (2022). PROTAC targeted protein degraders: the past is prologue. *Nat Rev Drug Discov*, 21, 181-200.

Jumper, J. et al. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature*, 596, 583-589.

Tsherniak, A. et al. (2017). Defining a cancer dependency map. *Cell*, 170, 564-576.

Tien, M. Z., Meyer, A. G., Sydykova, D. K., Spielman, S. J., & Wilke, C. O. (2013). Maximum allowed solvent accessibilities of residues in proteins. *PLoS One*, 8, e80635.

## Code and data availability

All code, data, and figures are available at https://github.com/crisprking/degradomap. The pipeline is reproducible end-to-end via the included CLI.
