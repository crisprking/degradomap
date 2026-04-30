"""
degradomap — empirical evaluation of public-data features for predicting
PROTAC E3 ligase tractability.

This package implements the pipeline described in the project's methods note:
  1. Query UniProt for human E3 ligases + verified PROTAC E3 positives
  2. Download AlphaFold structures and compute structural druggability features
  3. Pull DepMap CRISPR essentiality + expression summaries
  4. Train a leave-one-out cross-validated classifier and report AUC

The headline finding: structural + essentiality + expression features do NOT
discriminate the 12 published PROTAC E3 ligases from 814 other human E3
candidates (best LOO-CV AUC ~ 0.59). See docs/methods.md for full discussion.
"""

__version__ = "0.1.0"
