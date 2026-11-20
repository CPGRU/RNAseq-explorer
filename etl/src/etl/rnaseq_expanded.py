from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from etl.db_conn import get_pg_connection

# etl/src/etl/rnaseq_expanded_stats.py -> parents[3] is RNAseq-explorer/
BASE_DIR = Path(__file__).resolve().parents[3]
OUTPUT_DIR = BASE_DIR / "data" / "processed"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def fetch_expression_data() -> pd.DataFrame:
    conn = get_pg_connection()
    try:
        query = """
            SELECT
                g.sample_id,
                g.ensg_id,
                s.dex,
                s.celltype,
                g.expression_value::DOUBLE PRECISION AS expr
            FROM gene_expression g
            JOIN samples s ON g.sample_id = s.sample_id
            WHERE g.expression_value IS NOT NULL;
        """
        return pd.read_sql(query, conn)
    finally:
        conn.close()


def compute_distribution_stats(df: pd.DataFrame) -> pd.DataFrame:
    """3. Gene expression distribution, zero-inflation, and mean-variance metrics."""
    work = df.copy()
    work["log2_expr"] = np.log2(np.clip(work["expr"], a_min=0.0, a_max=None) + 1.0)
    work["is_zero"] = work["expr"] == 0.0

    grouped = work.groupby("ensg_id", as_index=False)

    dist_df = grouped.agg(
        n_samples=("expr", "count"),
        n_zero=("is_zero", "sum"),
        raw_mean=("expr", "mean"),
        raw_var=("expr", "var"),
        log2_mean=("log2_expr", "mean"),
        log2_var=("log2_expr", "var"),
        log2_q25=("log2_expr", lambda s: np.percentile(s, 25)),
        log2_q75=("log2_expr", lambda s: np.percentile(s, 75)),
    )

    dist_df["pct_zero"] = np.round((dist_df["n_zero"] / dist_df["n_samples"]) * 100.0, 2)
    dist_df["pct_nonzero"] = np.round(100.0 - dist_df["pct_zero"], 2)
    dist_df["log2_iqr"] = dist_df["log2_q75"] - dist_df["log2_q25"]

    return dist_df.drop(columns=["n_zero", "log2_q25", "log2_q75"])


def compute_paired_ttest_stats(df: pd.DataFrame) -> pd.DataFrame:
    """1. Paired t-test across matched cell types (treated vs. control)."""
    work = df.copy()
    work["log2_expr"] = np.log2(np.clip(work["expr"], a_min=0.0, a_max=None) + 1.0)

    # Pivot per gene and celltype to get control and treated side-by-side
    piv = (
        work.pivot_table(
            index=["ensg_id", "celltype"],
            columns="dex",
            values="log2_expr",
            aggfunc="mean",
        )
        .dropna(subset=["control", "treated"])
        .copy()
    )

    piv["paired_diff"] = piv["treated"] - piv["control"]

    rows = []
    for gene_id, grp in piv.groupby(level=0):
        n_pairs = len(grp)
        if n_pairs < 2:
            continue

        ctrl_vals = grp["control"].to_numpy(dtype=float)
        trt_vals = grp["treated"].to_numpy(dtype=float)
        diffs = trt_vals - ctrl_vals

        # Skip genes with identical differences across all pairs (zero variance -> NaN t-stat)
        if np.allclose(diffs, diffs[0]):
            continue

        t_stat, p_val = stats.ttest_rel(trt_vals, ctrl_vals)
        if not np.isfinite(t_stat) or not np.isfinite(p_val):
            continue

        rows.append({
            "ensg_id": gene_id,
            "n_pairs": int(n_pairs),
            "mean_log2_control": float(np.mean(ctrl_vals)),
            "mean_log2_treated": float(np.mean(trt_vals)),
            "paired_log2_fc": float(np.mean(diffs)),
            "paired_t_stat": float(t_stat),
            "paired_p_value": float(p_val),
        })

    paired_df = pd.DataFrame(rows)
    if paired_df.empty:
        return paired_df

    # Benjamini-Hochberg FDR adjustment across tested genes
    paired_df["paired_padj_bh"] = stats.false_discovery_control(
        paired_df["paired_p_value"].to_numpy(dtype=float),
        method="bh",
    )
    return paired_df


def add_volcano_and_ma_fields(merged_df: pd.DataFrame) -> pd.DataFrame:
    """4. Volcano plot and MA plot fields."""
    out = merged_df.copy()

    # Clamp padj away from 0 to avoid -log10(0) = inf
    safe_padj = np.clip(out["paired_padj_bh"].to_numpy(dtype=float), 1e-300, 1.0)
    out["neg_log10_padj"] = -np.log10(safe_padj)

    # MA plot coordinates: x = average abundance (log2_mean), y = fold change (paired_log2_fc)
    out["ma_x_mean_log2"] = out["log2_mean"]
    out["ma_y_log2_fc"] = out["paired_log2_fc"]

    # Volcano plot coordinates: x = fold change (paired_log2_fc), y = -log10(padj)
    out["volcano_x_log2_fc"] = out["paired_log2_fc"]
    out["volcano_y_neg_log10_padj"] = out["neg_log10_padj"]

    # Regulation classification (FDR < 0.05 and |log2FC| >= 1.0)
    is_up = (out["paired_padj_bh"] < 0.05) & (out["paired_log2_fc"] >= 1.0)
    is_down = (out["paired_padj_bh"] < 0.05) & (out["paired_log2_fc"] <= -1.0)

    out["regulation"] = np.where(is_up, "up", np.where(is_down, "down", "ns"))
    out["is_significant"] = is_up | is_down

    # Round numeric columns cleanly for readability and smaller JSON size
    float_cols = [
        "raw_mean",
        "raw_var",
        "log2_mean",
        "log2_var",
        "log2_iqr",
        "mean_log2_control",
        "mean_log2_treated",
        "paired_log2_fc",
        "paired_t_stat",
        "neg_log10_padj",
        "ma_x_mean_log2",
        "ma_y_log2_fc",
        "volcano_x_log2_fc",
        "volcano_y_neg_log10_padj",
    ]
    for col in float_cols:
        if col in out.columns:
            out[col] = np.round(out[col], 4)

    for p_col in ("paired_p_value", "paired_padj_bh"):
        if p_col in out.columns:
            out[p_col] = out[p_col].apply(lambda v: float(f"{v:.6e}"))

    return out.sort_values(by=["paired_padj_bh", "paired_p_value"], ascending=[True, True]).reset_index(
        drop=True
    )


def main() -> None:
    print("Fetching expression data from PostgreSQL...")
    raw_df = fetch_expression_data()
    print(f"Fetched {len(raw_df)} rows across {raw_df['ensg_id'].nunique()} genes.")

    dist_df = compute_distribution_stats(raw_df)
    paired_df = compute_paired_ttest_stats(raw_df)

    # Merge paired stats with distribution stats (inner join on genes with valid paired tests)
    combined_df = pd.merge(paired_df, dist_df, on="ensg_id", how="inner")
    final_df = add_volcano_and_ma_fields(combined_df)

    csv_path = OUTPUT_DIR / "rnaseq_paired_volcano_ma.csv"
    json_path = OUTPUT_DIR / "rnaseq_paired_volcano_ma.json"

    final_df.to_csv(csv_path, index=False)
    final_df.to_json(json_path, orient="records", indent=2)

    print(f"Saved CSV : {csv_path}")
    print(f"Saved JSON: {json_path}")
    print(f"Total genes in output: {len(final_df)}")
    print("Regulation counts:", final_df["regulation"].value_counts().to_dict())
    print("\nTop 5 genes preview:")
    preview_cols = [
        "ensg_id",
        "paired_log2_fc",
        "paired_p_value",
        "paired_padj_bh",
        "neg_log10_padj",
        "pct_zero",
        "regulation",
    ]
    print(final_df[preview_cols].head(5).to_string(index=False))


if __name__ == "__main__":
    main()