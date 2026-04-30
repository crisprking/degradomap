"""Command-line entry point: `python -m degradomap.cli`"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

from . import e3_list, structures, depmap, analysis


def main():
    ap = argparse.ArgumentParser(
        description="Run the degradomap pipeline end-to-end.")
    ap.add_argument("--workdir", type=Path, default=Path("./degradomap_run"))
    ap.add_argument("--skip-download", action="store_true",
                    help="Use existing files in workdir if present")
    args = ap.parse_args()
    W = args.workdir; W.mkdir(parents=True, exist_ok=True)

    print("=== Step 1: build E3 list ===")
    e3_path = W / "e3_ligases_verified.csv"
    if e3_path.exists() and args.skip_download:
        e3 = pd.read_csv(e3_path)
    else:
        e3 = e3_list.build_e3_list(out_path=e3_path)
    print(f"  {len(e3)} E3 candidates, {e3['protac_validated'].sum()} validated PROTAC E3s")

    print("=== Step 2: AlphaFold structures ===")
    pdb_dir = W / "pdb"
    structures.download_structures(e3["Entry"].tolist(), pdb_dir)
    n_pdb = len(list(pdb_dir.glob("*.pdb")))
    print(f"  {n_pdb} structures cached")

    print("=== Step 3: structural features ===")
    feat_path = W / "features.csv"
    rows = []
    for acc in e3["Entry"]:
        pdb = pdb_dir / f"{acc}.pdb"
        if not pdb.exists(): continue
        f = structures.compute_features(pdb)
        if f is None: continue
        f["accession"] = acc
        rows.append(f)
    feat_df = (pd.DataFrame(rows).merge(
        e3[["Entry","intended_gene","gene_symbol","protac_validated","Length"]],
        left_on="accession", right_on="Entry", how="left"))
    feat_df.to_csv(feat_path, index=False)
    print(f"  features for {len(feat_df)} proteins")

    print("=== Step 4: DepMap data ===")
    manifest = depmap.get_manifest()
    release = depmap.latest_release(manifest)
    print(f"  using DepMap release: {release}")
    ge_path = depmap.download_file(manifest, release, "CRISPRGeneEffect.csv", W / "gene_effect.csv")
    exp_path = depmap.download_file(manifest, release, "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", W / "expression_tpm.csv")

    print("=== Step 5: essentiality + expression summaries ===")
    e3_genes = set(feat_df["gene_symbol"].dropna())
    ess_df = depmap.essentiality_summary(ge_path, e3_genes)
    exp_df = depmap.expression_summary(exp_path, e3_genes)
    ess_df.to_csv(W / "e3_essentiality.csv", index=False)
    exp_df.to_csv(W / "e3_expression.csv", index=False)

    print("=== Step 6: merge and analyze ===")
    merged = (feat_df.merge(ess_df, on="gene_symbol", how="left")
                     .merge(exp_df, on="gene_symbol", how="left"))
    merged.to_csv(W / "merged_dataset.csv", index=False)

    result = analysis.run_experiment(merged)
    print(f"\nResult: n={result['n_total']} ({result['n_positives']} pos, {result['n_negatives']} neg)")
    print(f"AUC structural: {result['auc']['structural']:.3f}")
    print(f"AUC biological: {result['auc']['biological']:.3f}")
    print(f"AUC combined:   {result['auc']['combined']:.3f}")

    result["ranked"].to_csv(W / "ranked_predictions.csv", index=False)
    print(f"\nFull pipeline complete. Outputs in {W}")


if __name__ == "__main__":
    main()
