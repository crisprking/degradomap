# degradomap

**Empirical evaluation of public-data features for predicting PROTAC E3 ligase tractability.**

> **TL;DR.** We tested whether AlphaFold-derived structural druggability features,
> DepMap CRISPR essentiality, and cell-line expression breadth can discriminate
> the 12 published PROTAC E3 ligases from 814 other human E3 candidates.
> They cannot. Best leave-one-out cross-validated AUC was 0.587 (biological
> features alone); combining with structural features did not improve performance.
> Top-ranked predictions were dominated by core ubiquitin-pathway machinery
> (SKP1, ELOB, DDB1, E2 conjugating enzymes), not novel PROTAC E3 candidates.
>
> See [docs/methods.md](docs/methods.md) for the full methods note.

## What this is

A reproducible pipeline that:

1. Builds a curated list of human E3 ligase candidates from UniProt (~870 proteins),
   with 12–13 published PROTAC E3 ligases verified by gene-symbol lookup
2. Downloads AlphaFold-predicted structures for each
3. Computes 10 structural druggability features (compactness, hydrophobic surface
   patches, pocket-like residue clusters, AlphaFold pLDDT-weighted)
4. Pulls DepMap CRISPR gene-effect (Chronos) and expression (TPM) data for the
   current public release
5. Joins everything into a per-gene feature matrix and runs the discrimination
   experiment with leave-one-out CV

## What this is not

A working PROTAC E3 predictor. The result is a null finding, reported honestly.
This repository is most useful as:

- A reference implementation of the relevant data integrations (UniProt to AlphaFold
  to DepMap to unified per-gene feature matrix)
- A documented null result that other groups can build on without repeating the work
- A cautionary example of why structural features alone are insufficient for this
  question

## Install

Tested with Python 3.10 to 3.12. Requires ~1 GB of disk for AlphaFold structures
and DepMap matrices.

    git clone https://github.com/crisprking/degradomap.git
    cd degradomap
    pip install -e .

## Run the pipeline

    degradomap run --output-dir ./degradomap_run

This will (in order): query UniProt, download AlphaFold structures, compute
structural features, pull DepMap data, build the merged feature matrix, and
run the discrimination experiment. Total runtime ~30 minutes on a single core.

## Data files included

The data/ directory contains analysis-ready summary tables, all small (<1 MB):

- e3_ligases_verified.csv: curated E3 ligase candidates with PROTAC positives flagged
- features.csv: structural features per protein
- e3_essentiality.csv: DepMap CRISPR essentiality summaries
- e3_expression.csv: DepMap expression breadth summaries
- merged_dataset.csv: all features merged for the experiment

Raw AlphaFold structures and full DepMap matrices are not included (regenerable
from the pipeline; ~1.5 GB combined).

## Key result

See figures/experiment_result.png. The 12 PROTAC-validated E3s split into two
groups: about half rank in the top 300 (CRBN, MDM2, VHL, BIRC2, DCAF15, XIAP),
the other half rank below 500 (DCAF1, RNF4, RNF114, KEAP1, BIRC3, FEM1B). This
bimodality reflects the heterogeneity of the PROTAC E3 positive set: substrate
adapters versus direct E3 ligases, essential versus non-essential, broadly
versus narrowly expressed.

## Citation

Trueba, A. (2026). degradomap: empirical evaluation of public-data features
for predicting PROTAC E3 ligase tractability.
https://github.com/crisprking/degradomap

## License

MIT. See LICENSE.

## About the author

Abraham Trueba (UC Berkeley). I work on computational biology and drug discovery, with a focus on what public data can and cannot tell us about protein tractability for targeted protein degradation.

This work uses public data from UniProt, AlphaFold DB (EMBL-EBI), and DepMap
(Broad Institute). Cite their primary publications when using their data.
