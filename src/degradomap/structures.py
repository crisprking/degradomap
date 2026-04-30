"""AlphaFold structure download and structural feature extraction."""
from __future__ import annotations
import requests, time, numpy as np, pandas as pd
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from scipy.spatial.distance import cdist
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix

# Tien et al. 2013 empirical max side-chain SASA values (Å²)
MAX_SASA = {"ALA":129.0,"ARG":274.0,"ASN":195.0,"ASP":193.0,"CYS":167.0,
            "GLU":223.0,"GLN":225.0,"GLY":104.0,"HIS":224.0,"ILE":197.0,
            "LEU":201.0,"LYS":236.0,"MET":224.0,"PHE":240.0,"PRO":159.0,
            "SER":155.0,"THR":172.0,"TRP":285.0,"TYR":263.0,"VAL":174.0}
HYDROPHOBIC = {"ALA","VAL","LEU","ILE","MET","PHE","TRP","PRO","TYR"}
AROMATIC = {"PHE","TRP","TYR","HIS"}


def fetch_alphafold(uniprot_id: str) -> tuple[str | None, str]:
    """Fetch AlphaFold PDB via the prediction API. Returns (pdb_text, status)."""
    try:
        api = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
        meta = requests.get(api, timeout=20)
        if meta.status_code == 404: return None, "no_entry"
        meta.raise_for_status()
        entries = meta.json()
        if not entries: return None, "empty"
        pdb = requests.get(entries[0]["pdbUrl"], timeout=60)
        pdb.raise_for_status()
        return pdb.text, "ok"
    except Exception as e:
        return None, f"err_{type(e).__name__}"


def download_structures(accessions: list[str], pdb_dir: Path,
                         max_workers: int = 4) -> dict[str, str]:
    """Download AlphaFold structures for a list of accessions, with caching."""
    pdb_dir = Path(pdb_dir); pdb_dir.mkdir(parents=True, exist_ok=True)

    def one(acc):
        p = pdb_dir / f"{acc}.pdb"
        if p.exists() and p.stat().st_size > 1000:
            return acc, "cached"
        text, status = fetch_alphafold(acc)
        if text: p.write_text(text)
        return acc, status

    results = {}
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        for fut in as_completed({ex.submit(one, a): a for a in accessions}):
            acc, st = fut.result()
            results[acc] = st
    return results


def compute_features(pdb_path: Path) -> dict | None:
    """Compute structural druggability features from a single PDB file.

    Returns a dict of features, or None if the protein is too small for analysis.
    """
    try:
        s = PDBParser(QUIET=True).get_structure("X", str(pdb_path))
        ShrakeRupley().compute(s, level="R")
        rows = []
        for ch in s[0]:
            for res in ch:
                if res.id[0] != " ": continue
                rn = res.get_resname()
                if rn not in MAX_SASA or "CA" not in res: continue
                rsasa = res.sasa / MAX_SASA[rn]
                rows.append({
                    "resname": rn, "rel_sasa": rsasa,
                    "plddt": float(np.mean([a.bfactor for a in res])),
                    "ca": res["CA"].coord,
                    "is_hp": rn in HYDROPHOBIC,
                    "is_arom": rn in AROMATIC,
                    "is_surface": rsasa > 0.20,
                    "is_pocket_like": 0.05 < rsasa < 0.40,
                })
        if len(rows) < 30: return None
        df = pd.DataFrame(rows)
        conf = df[df["plddt"] > 70].reset_index(drop=True)
        if len(conf) < 30: return None
        n = len(conf)

        hp_idx = conf[conf["is_hp"] & conf["is_surface"]].index
        max_patch, n_patches5 = 0, 0
        if len(hp_idx) >= 2:
            coords = np.vstack(conf.loc[hp_idx, "ca"].values)
            adj = (cdist(coords, coords) <= 7.0).astype(int)
            np.fill_diagonal(adj, 0)
            _, labels = connected_components(csr_matrix(adj), directed=False)
            sizes = pd.Series(labels).value_counts()
            max_patch = int(sizes.max())
            n_patches5 = int((sizes >= 5).sum())

        pl = conf[conf["is_hp"] & conf["is_pocket_like"]]
        max_pocket = 0
        if len(pl) >= 3:
            pc = np.vstack(pl["ca"].values)
            adj = (cdist(pc, pc) <= 8.0).astype(int)
            np.fill_diagonal(adj, 0)
            _, labels = connected_components(csr_matrix(adj), directed=False)
            max_pocket = int(pd.Series(labels).value_counts().max())

        coords_all = np.vstack(conf["ca"].values)
        rg = float(np.sqrt(((coords_all - coords_all.mean(0)) ** 2).sum(1).mean()))
        surf = conf[conf["is_surface"]]

        return {
            "n_residues_total": len(df),
            "n_residues_confident": n,
            "confident_fraction": n / len(df),
            "mean_plddt": float(conf["plddt"].mean()),
            "compactness": rg / np.sqrt(n),
            "hp_surf_frac": float(surf["is_hp"].mean()) if len(surf) else 0.0,
            "arom_surf_frac": float(surf["is_arom"].mean()) if len(surf) else 0.0,
            "max_hp_patch": max_patch,
            "n_hp_patches_5plus": n_patches5,
            "pocket_density": len(pl) / n,
            "max_pocket_cluster": max_pocket,
            "pocket_cluster_per_100res": max_pocket / n * 100,
        }
    except Exception:
        return None
