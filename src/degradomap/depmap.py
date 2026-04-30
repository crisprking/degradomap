"""DepMap data download and per-gene essentiality/expression summaries."""
from __future__ import annotations
import requests, re, pandas as pd
from io import StringIO
from pathlib import Path

MANIFEST_URL = "https://depmap.org/portal/api/download/files"
HEADERS = {"User-Agent": "Mozilla/5.0 degradomap/0.1"}


def get_manifest() -> pd.DataFrame:
    """Fetch the DepMap file manifest as a DataFrame."""
    r = requests.get(MANIFEST_URL, headers=HEADERS, timeout=60)
    r.raise_for_status()
    return pd.read_csv(StringIO(r.text))


def latest_release(manifest: pd.DataFrame) -> str:
    """Return the most recent release name in the manifest."""
    manifest = manifest.copy()
    manifest["release_date"] = pd.to_datetime(manifest["release_date"])
    return manifest.sort_values("release_date", ascending=False)["release"].iloc[0]


def download_file(manifest: pd.DataFrame, release: str, filename: str,
                  out_path: Path) -> Path:
    """Download a specific file from a specific release."""
    out_path = Path(out_path)
    if out_path.exists() and out_path.stat().st_size > 1000:
        return out_path
    row = manifest[(manifest["release"] == release)
                   & (manifest["filename"] == filename)]
    if row.empty:
        raise ValueError(f"{filename} not found in release {release}")
    url = row["url"].iloc[0]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, headers=HEADERS, stream=True, timeout=600) as resp:
        resp.raise_for_status()
        with open(out_path, "wb") as f:
            for chunk in resp.iter_content(chunk_size=1024 * 1024):
                f.write(chunk)
    return out_path


_COL_PATTERN = re.compile(r"^([A-Z0-9_-]+)\s*\(\d+\)\s*$")


def parse_gene_columns(columns) -> dict[str, str]:
    """Map DepMap 'GENE (entrez_id)' headers to bare gene symbols."""
    out = {}
    for c in columns:
        m = _COL_PATTERN.match(c)
        if m:
            out[c] = m.group(1)
    return out


def essentiality_summary(gene_effect_csv: Path, gene_symbols: set[str]) -> pd.DataFrame:
    """Compute per-gene essentiality stats for a set of gene symbols.

    Returns DataFrame with columns:
      gene_symbol, n_lines_screened, mean_chronos, median_chronos,
      frac_essential (Chronos < -0.5), frac_strongly_essential (Chronos < -1.0)
    """
    header = pd.read_csv(gene_effect_csv, nrows=0)
    col_map = parse_gene_columns(header.columns)
    sym_to_col = {sym: col for col, sym in col_map.items()}
    matched = gene_symbols & set(sym_to_col)
    cols = [header.columns[0]] + [sym_to_col[g] for g in matched]
    df = pd.read_csv(gene_effect_csv, usecols=cols)
    df = df.rename(columns={**{sym_to_col[g]: g for g in matched},
                            df.columns[0]: "DepMap_ID"})
    rows = []
    for g in matched:
        s = df[g].dropna()
        if len(s) < 50: continue
        rows.append({
            "gene_symbol": g,
            "n_lines_screened": len(s),
            "mean_chronos": float(s.mean()),
            "median_chronos": float(s.median()),
            "frac_essential": float((s < -0.5).mean()),
            "frac_strongly_essential": float((s < -1.0).mean()),
        })
    return pd.DataFrame(rows)


def expression_summary(expression_csv: Path, gene_symbols: set[str]) -> pd.DataFrame:
    """Compute per-gene expression breadth stats.

    Input is the DepMap log2(TPM+1) matrix.
    Returns DataFrame with columns:
      gene_symbol, n_lines_expression, mean_log_tpm, median_log_tpm,
      frac_expressed (log2(TPM+1) > 1, i.e. TPM > 1),
      frac_high_expressed (log2(TPM+1) > 3, i.e. TPM > 7)
    """
    header = pd.read_csv(expression_csv, nrows=0)
    col_map = parse_gene_columns(header.columns)
    sym_to_col = {sym: col for col, sym in col_map.items()}
    matched = gene_symbols & set(sym_to_col)
    cols = [header.columns[0]] + [sym_to_col[g] for g in matched]
    df = pd.read_csv(expression_csv, usecols=cols)
    df = df.rename(columns={**{sym_to_col[g]: g for g in matched},
                            df.columns[0]: "DepMap_ID"})
    rows = []
    for g in matched:
        s = df[g].dropna()
        if len(s) < 50: continue
        rows.append({
            "gene_symbol": g,
            "n_lines_expression": len(s),
            "mean_log_tpm": float(s.mean()),
            "median_log_tpm": float(s.median()),
            "frac_expressed": float((s > 1.0).mean()),
            "frac_high_expressed": float((s > 3.0).mean()),
        })
    return pd.DataFrame(rows)
