# degradomap

**Empirical evaluation of public-data features for predicting PROTAC E3 ligase tractability — a living null result.**

Five pre-registered experiments. Five failures. Documented honestly.

> **TL;DR (v5).** Biological features (DepMap essentiality + expression) dominate our
> classifier and create a UPS-pathway shortcut: the model learns "essential + broadly expressed"
> rather than "covalently ligandable." Structural features alone perform barely better than chance
> (LOO AUC ~0.55). None of our four pre-registered v5 tests passed. The covalent PROTAC E3
> problem remains unsolved by public-data features.
>
> See [docs/methods.md](docs/methods.md) for the full methods note and version-by-version results.

## Version history

| Version | Question | Pre-registered | Result |
|---------|----------|----------------|--------|
| v1 | Can AlphaFold structural features + DepMap data discriminate PROTAC E3s from the human E3 universe? | No (retrospective) | Null. Best LOO AUC 0.587 (bio features only). |
| v2 | Does mechanism-stratified analysis (pocket/glue/covalent) rescue signal? | Yes | Null. All three sub-experiments fail independently. |
| v3 | Does held-out external validation on 5 post-2024 E3s confirm null? | Yes | Confirmed null. 3/5 held-out E3s rank below median. |
| v4 | Do propeller-fold features (DCAF11 structural homologs) add signal? | Yes | Catastrophic failure. WD40 propeller fold shortcut. |
| v5 | With a diversified covalent positive set (6 proteins, 4 folds), does anything work? | Yes | Null + UPS-pathway shortcut diagnosed. 1/4 tests passed. |

## What this repo is

A reproducible, version-controlled record of a systematic null result in computational target identification for targeted protein degradation. Includes:

- Full pipeline (UniProt → AlphaFold → DepMap → feature matrix → LOO classifier)
- Pre-registration documents for v2–v5 ([docs/preregistrations/](docs/preregistrations/))
- Version-by-version data files in [data/](data/)
- Source modules for each analytical approach, including v4 (documented failure)

## What this repo is not

A working PROTAC E3 predictor. This is the honest baseline that shows why naive public-data approaches fail and what the next step actually requires.

## Install

Tested with Python 3.10–3.12.

```bash
git clone https://github.com/crisprking/degradomap.git
cd degradomap
pip install -e .
```

## Run the pipeline

```bash
degradomap run --output-dir ./degradomap_run
```

Queries UniProt, downloads AlphaFold structures, computes features, pulls DepMap data, runs the discrimination experiment. ~30 min on a single core.

## Data files

| File | Description |
|------|-------------|
| `data/e3_ligases_verified.csv` | Human E3 candidates with PROTAC positives flagged |
| `data/features.csv` | Structural druggability features per protein |
| `data/merged_dataset.csv` | v1 unified feature matrix |
| `data/per_mechanism_summary.csv` | v2 mechanism-split results |
| `data/external_validation_results.csv` | v3 held-out E3 scores |
| `data/propeller_features.csv` | v4 propeller fold features (failed approach) |
| `data/v5_loo_coherence_results.csv` | v5 LOO by covalent positive |
| `data/v5_ablation_loo.csv` | v5 full vs structural-only ablation |
| `data/v5_structural_top15_audit.csv` | v5 manual audit of structural top-15 |

## Key failure modes diagnosed

1. **UPS-pathway shortcut** — DepMap essentiality/expression features correlate with UPS membership, not covalent ligandability. The classifier learns the shortcut.
2. **Fold-family shortcut** — propeller structural features (v4) caused the model to over-rank all WD40 propeller proteins regardless of covalent site chemistry.
3. **Positive set sparsity** — 6 verified covalent positives out of ~870 candidates is too sparse for tabular feature generalization.

## Citation

Trueba, A. (2026). degradomap: empirical evaluation of public-data features for predicting PROTAC E3 ligase tractability. https://github.com/crisprking/degradomap

## License

MIT. See [LICENSE](LICENSE).

## About the author

Abraham Trueba (UC Berkeley). Computational biology and drug discovery, focused on what public data can and cannot tell us about protein tractability for targeted protein degradation.

Uses public data from UniProt, AlphaFold DB (EMBL-EBI), and DepMap (Broad Institute). Cite their primary publications when using their data.
