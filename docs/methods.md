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

---

## §4 — v2: Mechanism-stratified analysis

**Pre-registered:** Yes (see docs/preregistrations/v2.md)

### Motivation

The v1 positive set is mechanistically heterogeneous: pocket binders (VHL, MDM2, KEAP1, XIAP/BIRC2/BIRC3), molecular glues (CRBN, DCAF15, DCAF16), and covalent binders (RNF4, RNF114, FEM1B, DCAF1). Pooling these into one positive class erases structure-function relationships.

### Classification

| Mechanism | Members |
|-----------|---------|
| pocket_binder | VHL, MDM2, KEAP1, XIAP, BIRC2, BIRC3 |
| molecular_glue | CRBN, DCAF15, DCAF16 |
| covalent | RNF4, RNF114, FEM1B, DCAF1 |

### Results

All three mechanism-specific experiments failed to produce better-than-chance discrimination. The pocket-binder sub-experiment showed the weakest signal (LOO AUC ~0.54); the covalent sub-experiment showed the strongest but still non-significant result (LOO AUC ~0.61). Permutation testing confirmed all results are consistent with chance.

---

## §5 — v3: External held-out validation

**Pre-registered:** Yes (see docs/preregistrations/v3.md)

Five PROTAC E3s published after data collection were held out and scored by the v1 model without retraining:

| Gene | Year | Mechanism | Rank | Pass? |
|------|------|-----------|------|-------|
| DCAF11 | 2024 | covalent | ~350 | No |
| RNF126 | 2024 | covalent | ~280 | Partial |
| FBXO22 | 2024 | covalent | ~420 | No |
| TRIM21 | 2023 | pocket | ~190 | Yes |
| STUB1 | 2023 | pocket | ~210 | Yes |

Pocket binders generalize weakly; covalent binders do not. Confirms v1 null with mechanism-specific detail.

---

## §6 — v4: Propeller-fold features (failed approach)

**Pre-registered:** Yes (see docs/preregistrations/v4.md)

### Motivation

DCAF11, a newly validated covalent PROTAC E3, is a WD40 β-propeller. Our existing structural features (compactness, hydrophobic patches, pocket clusters) are not informative for propeller folds, which lack deep druggable pockets but present shallow surface pockets on the top face.

### Approach

Added 4 propeller-fold-specific features: propeller_score (HHsearch against β-propeller SCOP family), top_face_hydrophobicity, top_face_cleft_depth, propeller_cys_proximity. All computed from AlphaFold structures.

### Result

**Catastrophic failure.** The propeller features caused the model to rank all WD40 propeller proteins highly (~800 proteins), generating massive false positives. DCAF11 LOO rank improved from ~350 to ~85, but SKP1, WD45, RACK1, and other non-E3 propeller proteins now rank in the top 20.

### Lesson

Fold-specific features create fold-family shortcuts. A model that learns "propeller-like → high score" cannot discriminate covalently-ligandable propeller E3s from the broad WD40 family. Fold-specificity requires residue-level cysteine reactivity modeling, not family-level structural features.

Code is preserved in `src/degradomap/propeller.py` for documentation.

---

## §7 — v5: Diversified covalent positive set + shortcut diagnosis

**Pre-registered:** Yes (see docs/preregistrations/v5.md)

### Motivation

The v2 covalent training set (n=4) was structurally homogeneous: 3 small RING-type proteins (RNF4, RNF114, FEM1B) plus DCAF1. This homogeneity allowed v4 to find the fold-family shortcut. v5 expands to 6 verified covalent positives across 4 folds.

### Positive set

| Gene | Accession | Fold | Citation |
|------|-----------|------|----------|
| RNF4 | P78317 | RING | Ward et al. Cell Chem Biol 2019 |
| RNF114 | Q9Y508 | RING | Spradlin et al. Nat Chem Biol 2019 |
| DCAF1 | Q9Y4B6 | multi-domain | Multiple 2022–2024 |
| DCAF11 | Q8TEB1 | WD40 propeller | Tin et al. BMCL 2024 |
| RNF126 | Q9BV68 | RING | Lim et al. ACS Cent Sci 2024 |
| FBXO22 | Q8NEZ5 | F-box/LRR | Nie et al. Nat Chem Biol 2024 |

FEM1B (Q92545) and DCAF16 (Q9NXF7) excluded: missing DepMap features and structural disorder, respectively. All accessions verified by gene-symbol lookup at runtime.

### Pre-registered tests

| # | Test | Threshold | Result |
|---|------|-----------|--------|
| 1 (PRIMARY) | Median LOO percentile ≤ 25% | Median ≤ 25 | **FAIL** — median 30.4% |
| 2 (SECONDARY) | Multi-domain Spearman ρ > 0.15, p < 0.05 | Both conditions | **FAIL** — ρ = 0.11, p = 0.18 |
| 3 (TERTIARY) | DCAF1 structural-only LOO ≤ 50% | Percentile ≤ 50 | **PARTIAL PASS** — percentile = 48 |
| 4 (QUATERNARY) | Structural top-15 contains ≥ 10 confirmed E3s | Count ≥ 10 | **FAIL** — 7 confirmed E3s, 5 UPS scaffold proteins |

### Shortcut diagnosis

Ablation experiment (full model vs structural-only model):

| Gene | Full model %ile | Structural-only %ile | Δ |
|------|----------------|---------------------|---|
| RNF4 | 8 | 34 | -26 |
| RNF114 | 12 | 41 | -29 |
| DCAF1 | 48 | 48 | 0 |
| DCAF11 | 71 | 58 | +13 |
| RNF126 | 65 | 55 | +10 |
| FBXO22 | 82 | 69 | +13 |

Biological features help RING-type positives (which are essential and broadly expressed) but hurt newer covalent positives (which are non-essential and narrowly expressed). This is the UPS-pathway shortcut: the model learns "essential + broadly expressed = UPS core member = positive" rather than "structurally covalent-reactive = positive."

### Structural top-15 audit

Manual review of top-15 structural-only predictions: 7 confirmed E3 ligases (in addition to training positives), 5 UPS scaffold proteins (SKP1, ELOB, DDB1, CUL5, RBBP7), 3 non-ubiquitin proteins (HDAC8, SETD7, PRMT5). The scaffold proteins appear due to AlphaFold structural similarity to RING-type folds, not because of covalent ligandability.

### Conclusion

Public-data features (structural + DepMap) cannot solve the covalent PROTAC E3 discrimination problem. The next step requires warhead-side data: CysDB cysteine reactivity scores, competitive ABPP datasets, or AlphaFold-predicted cysteine accessibility + electrophilicity proxies.
