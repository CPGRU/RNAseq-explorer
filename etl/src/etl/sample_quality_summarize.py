from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from etl.db_conn import get_pg_connection

BASE_DIR = Path(__file__).resolve().parents[3]
OUTPUT_DIR = BASE_DIR / "data" / "processed"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def fetch_data() -> pd.DataFrame:
    conn = get_pg_connection()
    try:
        query = """
            SELECT
                g.sample_id,
                g.ensg_id,
                s.dex,
                s.celltype,
                s.geo_id,
                g.expression_value::DOUBLE PRECISION AS expr
            FROM gene_expression g
            JOIN samples s ON g.sample_id = s.sample_id
            WHERE g.expression_value IS NOT NULL;
        """
        return pd.read_sql(query, conn)
    finally:
        conn.close()


def compute_sample_pca_and_corr(df: pd.DataFrame) -> tuple[list[dict], dict]:
    df = df.copy()
    df["log2_expr"] = np.log2(df["expr"] + 1.0)

    # Wide matrix: genes x samples
    wide = df.pivot(index="ensg_id", columns="sample_id", values="log2_expr").dropna()

    # 1. Sample correlation matrix (samples x samples)
    corr_df = wide.corr(method="pearson").round(4)
    corr_payload = {
        "samples": list(corr_df.columns),
        "matrix": corr_df.values.tolist(),
    }

    # 2. PCA via SVD on centered log2 expression (samples as rows -> transpose)
    X = wide.T.values  # shape: (n_samples, n_genes)
    X_centered = X - X.mean(axis=0, keepdims=True)
    U, S, _ = np.linalg.svd(X_centered, full_matrices=False)
    coords = U[:, :2] * S[:2]
    explained_var = (S**2) / np.sum(S**2)

    meta = (
        df[["sample_id", "dex", "celltype", "geo_id"]]
        .drop_duplicates("sample_id")
        .set_index("sample_id")
        .loc[wide.columns]
        .reset_index()
    )

    pca_records = []
    for idx, row in meta.iterrows():
        pca_records.append({
            "sample_id": row["sample_id"],
            "dex": row["dex"],
            "celltype": row["celltype"],
            "geo_id": row["geo_id"],
            "pc1": round(float(coords[idx, 0]), 4),
            "pc2": round(float(coords[idx, 1]), 4),
            "pc1_var_explained": round(float(explained_var[0]), 4),
            "pc2_var_explained": round(float(explained_var[1]), 4),
        })

    return pca_records, corr_payload


def compute_paired_gene_stats(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["log2_expr"] = np.log2(df["expr"] + 1.0)

    # Pivot to (ensg_id, celltype) with columns 'control' and 'treated'
    piv = df.pivot_table(
        index=["ensg_id", "celltype"],
        columns="dex",
        values="log2_expr",
        aggfunc="mean",
    ).dropna(subset=["control", "treated"])

    piv["diff"] = piv["treated"] - piv["control"]

    rows = []
    for gene_id, grp in piv.groupby(level=0):
        # Requires paired observations across cell types
        if len(grp) < 2:
            continue
        diffs = grp["diff"].values
        if np.allclose(diffs, diffs[0]):
            continue

        t_stat, p_val = stats.ttest_rel(grp["treated"].values, grp["control"].values)
        mean_diff = float(np.mean(diffs))
        rows.append({
            "ensg_id": gene_id,
            "n_pairs": int(len(grp)),
            "paired_mean_log2_diff": round(mean_diff, 4),
            "paired_t_stat": round(float(t_stat), 4) if np.isfinite(t_stat) else None,
            "paired_p_value": float(p_val) if np.isfinite(p_val) else None,
        })

    out = pd.DataFrame(rows).dropna(subset=["paired_p_value"]).copy()
    out["paired_padj_bh"] = stats.false_discovery_control(out["paired_p_value"], method="bh")
    out["neg_log10_padj"] = -np.log10(np.clip(out["paired_padj_bh"], 1e-300, 1.0))
    out["regulation"] = np.where(
        (out["paired_padj_bh"] < 0.05) & (out["paired_mean_log2_diff"] >= 1.0),
        "up",
        np.where(
            (out["paired_padj_bh"] < 0.05) & (out["paired_mean_log2_diff"] <= -1.0),
            "down",
            "ns",
        ),
    )
    return out.sort_values("paired_padj_bh").reset_index(drop=True)


def main() -> None:
    df = fetch_data()
    pca_records, corr_payload = compute_sample_pca_and_corr(df)
    paired_df = compute_paired_gene_stats(df)

    (OUTPUT_DIR / "sample_pca.json").write_text(
        pd.Series(pca_records).to_json(orient="values", indent=2)
    )
    import json
    (OUTPUT_DIR / "sample_corr.json").write_text(json.dumps(corr_payload, indent=2))
    paired_df.to_json(OUTPUT_DIR / "paired_gene_stats.json", orient="records", indent=2)

    summary = {
        "total_genes_tested": int(len(paired_df)),
        "upregulated": int((paired_df["regulation"] == "up").sum()),
        "downregulated": int((paired_df["regulation"] == "down").sum()),
        "non_significant": int((paired_df["regulation"] == "ns").sum()),
    }
    (OUTPUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2))

    print("Exported:")
    print("  - data/processed/summary.json")
    print("  - data/processed/sample_pca.json")
    print("  - data/processed/sample_corr.json")
    print("  - data/processed/paired_gene_stats.json")
    print("Summary:", summary)


if __name__ == "__main__":
    main()