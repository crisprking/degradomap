"""Build a curated human E3 ligase list with verified PROTAC positives."""
from __future__ import annotations
import requests, time, pandas as pd
from io import StringIO
from pathlib import Path

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb/search"

# Validated PROTAC E3 ligases — gene symbols. Resolved to UniProt accessions
# at runtime via gene-symbol search to avoid stale hardcoded accessions.
PROTAC_E3_GENES = [
    "CRBN", "VHL", "XIAP", "BIRC2", "BIRC3", "MDM2",
    "DCAF15", "DCAF16", "DCAF1", "RNF114", "RNF4", "KEAP1", "FEM1B",
]

E3_QUERIES = {
    "direct_e3": (
        "(organism_id:9606) AND (reviewed:true) AND "
        "(cc_function:\"E3 ubiquitin\" OR cc_function:\"ubiquitin ligase\" "
        "OR go:0061630 OR go:0004842)"
    ),
    "substrate_adapters": (
        "(organism_id:9606) AND (reviewed:true) AND "
        "(cc_function:\"substrate receptor\" OR cc_function:\"DDB1\" "
        "OR cc_function:\"CUL4\" OR protein_name:DCAF OR go:1990756)"
    ),
    "f_box": (
        "(organism_id:9606) AND (reviewed:true) AND "
        "(protein_name:\"F-box\" OR cc_similarity:\"F-box family\")"
    ),
    "btb_socs": (
        "(organism_id:9606) AND (reviewed:true) AND "
        "(protein_name:\"SOCS box\" OR protein_name:\"BTB/POZ\" "
        "OR cc_function:\"Cullin-3\")"
    ),
}

FIELDS = "accession,id,gene_names,protein_name,length"


def _paginated_query(query: str) -> pd.DataFrame:
    """Run a paginated UniProt query, return concatenated DataFrame."""
    url = UNIPROT_BASE
    params = {"query": query, "fields": FIELDS, "format": "tsv", "size": 500}
    dfs, page = [], 0
    while True:
        page += 1
        r = (requests.get(url, params=params, timeout=60) if page == 1
             else requests.get(url, timeout=60))
        r.raise_for_status()
        df = pd.read_csv(StringIO(r.text), sep="\t")
        if len(df) == 0: break
        dfs.append(df)
        if "next" in r.links:
            url = r.links["next"]["url"]
            time.sleep(0.4)
        else:
            break
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def _resolve_gene(symbol: str) -> dict | None:
    """Look up a gene symbol's reviewed Swiss-Prot entry."""
    params = {
        "query": f"(gene_exact:{symbol}) AND (organism_id:9606) AND (reviewed:true)",
        "fields": FIELDS, "format": "tsv", "size": 5,
    }
    r = requests.get(UNIPROT_BASE, params=params, timeout=30)
    r.raise_for_status()
    df = pd.read_csv(StringIO(r.text), sep="\t")
    if len(df) == 0: return None
    df["primary"] = df["Gene Names"].astype(str).str.split().str[0]
    match = df[df["primary"] == symbol]
    return (match.iloc[0] if len(match) else df.iloc[0]).to_dict()


def build_e3_list(out_path: Path | str | None = None) -> pd.DataFrame:
    """Build the full E3 ligase candidate list with verified PROTAC positives.

    Returns a DataFrame with columns: Entry, Entry Name, Gene Names,
    Protein names, Length, source_query, intended_gene, protac_validated,
    gene_symbol.
    """
    parts = []
    for label, q in E3_QUERIES.items():
        df = _paginated_query(q)
        df["source_query"] = label
        parts.append(df)
    e3_raw = pd.concat(parts, ignore_index=True)

    e3 = (e3_raw.groupby("Entry")
          .agg({**{c: "first" for c in e3_raw.columns
                   if c not in ["Entry", "source_query"]},
                "source_query": lambda x: ",".join(sorted(set(x)))})
          .reset_index())

    # Resolve PROTAC E3 positives by gene symbol
    extras = []
    for g in PROTAC_E3_GENES:
        info = _resolve_gene(g)
        if info is None: continue
        extras.append({
            **{k: info.get(k, "") for k in
               ["Entry", "Entry Name", "Gene Names", "Protein names", "Length"]},
            "source_query": "manual_protac_e3",
            "intended_gene": g,
        })
        time.sleep(0.2)
    extras_df = pd.DataFrame(extras)

    e3 = pd.concat([e3.assign(intended_gene=None), extras_df], ignore_index=True)
    e3["protac_validated"] = e3["Entry"].isin(extras_df["Entry"])
    e3 = (e3.sort_values("protac_validated", ascending=False)
          .drop_duplicates("Entry").reset_index(drop=True))

    protac_map = dict(zip(extras_df["Entry"], extras_df["intended_gene"]))
    e3.loc[e3["Entry"].isin(protac_map), "intended_gene"] = (
        e3.loc[e3["Entry"].isin(protac_map), "Entry"].map(protac_map))

    e3["gene_symbol"] = e3["Gene Names"].astype(str).str.split().str[0]

    if out_path is not None:
        e3.to_csv(out_path, index=False)
    return e3
