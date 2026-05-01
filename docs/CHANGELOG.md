# Changelog

## v5 — Diversified covalent positive set + shortcut diagnosis
*May 2026*

- Expanded covalent positive set from 4 to 6 proteins across 4 structural folds
- Added `src/degradomap/covalent.py`: LOO coherence, ablation, and universe scoring for v5
- Pre-registered 4 tests before running analysis; 1 partial pass, 3 fails
- Diagnosed UPS-pathway shortcut driven by DepMap essentiality/expression features
- Manual audit of structural top-15 predictions
- Data: `data/v5_loo_coherence_results.csv`, `v5_ablation_loo.csv`, `v5_top20_full.csv`, `v5_top20_structural.csv`, `v5_structural_top15_audit.csv`

## v4 — Propeller-fold features (failed)
*April 2026*

- Added propeller-fold-specific structural features targeting DCAF11 (WD40 covalent E3)
- Pre-registered before running; result was catastrophic (fold-family shortcut)
- Code preserved in `src/degradomap/propeller.py` for documentation
- Data: `data/propeller_features.csv`, `data/merged_dataset_v4.csv`

## v3 — External held-out validation
*April 2026*

- Held out 5 post-2024 PROTAC E3s and scored with v1 model without retraining
- Pre-registered; confirmed v1 null with mechanism-specific detail
- Pocket binders generalize weakly; covalent binders do not
- Data: `data/external_validation_results.csv`

## v2 — Mechanism-stratified analysis
*April 2026*

- Split 12 PROTAC E3 positives into pocket_binder, molecular_glue, covalent classes
- Pre-registered; all three sub-experiments fail independently
- Added `src/degradomap/mechanism.py`
- Data: `data/per_mechanism_summary.csv`, `data/per_mechanism_feature_tests.csv`

## v1 — Initial analysis
*April 2026*

- Built pipeline: UniProt → AlphaFold → DepMap → feature matrix → LOO classifier
- Best LOO AUC 0.587 (biological features alone); structural features do not improve
- Top-ranked predictions dominated by UPS core machinery (SKP1, ELOB, DDB1)
- Retrospective; no pre-registration
