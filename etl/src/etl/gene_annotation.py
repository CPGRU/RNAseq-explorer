import pandas as pd
import requests


def fetch_ensembl_annotations(
    ensg_ids: list[str], batch_size: int = 1000
) -> dict[str, dict[str, str | None]]:
    """Fetch gene symbol and full gene name for clean Ensembl IDs (e.g., ENSG00000084710)."""
    clean_ids = sorted(
        {str(gid).strip() for gid in ensg_ids if pd.notna(gid) and str(gid).strip()}
    )
    annot_map: dict[str, dict[str, str | None]] = {}

    url = "https://mygene.info/v3/gene"
    for start in range(0, len(clean_ids), batch_size):
        chunk = clean_ids[start : start + batch_size]
        resp = requests.post(
            url,
            json={"ids": chunk, "fields": "symbol,name"},
            timeout=30,
        )
        resp.raise_for_status()
        records = resp.json()
        if isinstance(records, dict):
            records = [records]

        for rec in records:
            if not isinstance(rec, dict) or rec.get("notfound"):
                continue
            qid = str(rec.get("query", "")).strip()
            if not qid:
                continue
            # Keep first non-empty symbol/name if multiple records share a query ID
            entry = annot_map.setdefault(qid, {"gene_symbol": None, "gene_name": None})
            if entry["gene_symbol"] is None and rec.get("symbol"):
                entry["gene_symbol"] = str(rec["symbol"])
            if entry["gene_name"] is None and rec.get("name"):
                entry["gene_name"] = str(rec["name"])

    return annot_map


def annotate_ensembl_df(
    df: pd.DataFrame, id_col: str = "ensg_id"
) -> pd.DataFrame:
    """Insert `gene_symbol` and `gene_name` columns right after `id_col`."""
    out = df.copy()
    ids = out[id_col].astype(str).str.strip()
    annot_map = fetch_ensembl_annotations(ids.tolist())

    symbols = ids.map(lambda x: annot_map.get(x, {}).get("gene_symbol"))
    names = ids.map(lambda x: annot_map.get(x, {}).get("gene_name"))

    insert_at = list(out.columns).index(id_col) + 1
    out.insert(insert_at, "gene_symbol", symbols)
    out.insert(insert_at + 1, "gene_name", names)
    return out