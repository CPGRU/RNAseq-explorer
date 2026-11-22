from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from etl.db_conn import get_pg_connection

# Project root: etl/src/etl/dge_analysis.py -> 3 levels up is RNAseq-explorer/
BASE_DIR = Path(__file__).resolve().parents[3]
OUTPUT_DIR = BASE_DIR / "data" / "processed"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def fetch_joined_expression() -> pd.DataFrame:
    conn = get_pg_connection()
    try:
        query = """
            SELECT
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


def run_analysis(df: pd.DataFrame) -> pd.DataFrame:
    # Pseudocount (+1) log2 transformation for exploratory comparison
    df = df.copy()
    df["log2_expr"] = np.log2(df["expr"].astype(float) + 1.0)

    rows = []
    for gene_id, grp in df.groupby("ensg_id"):
        treated_raw = grp.loc[grp["dex"] == "treated", "expr"]
        control_raw = grp.loc[grp["dex"] == "control", "expr"]

        treated_log2 = grp.loc[grp["dex"] == "treated", "log2_expr"]
        control_log2 = grp.loc[grp["dex"] == "control", "log2_expr"]

        if len(treated_log2) < 2 or len(control_log2) < 2:
            continue

        mean_control = float(control_raw.mean())
        mean_treated = float(treated_raw.mean())
        mean_log2_control = float(control_log2.mean())
        mean_log2_treated = float(treated_log2.mean())

        # Log2 fold change between treated and control log2-transformed means
        log2_fc = mean_log2_treated - mean_log2_control

        # Welch's two-sample t-test
        _, p_val = stats.ttest_ind(treated_log2, control_log2, equal_var=False)

        rows.append({
            "ensg_id": gene_id,
            "mean_control": mean_control,
            "mean_treated": mean_treated,
            "mean_log2_control": mean_log2_control,
            "mean_log2_treated": mean_log2_treated,
            "log2_fold_change": float(log2_fc),
            "p_value": float(p_val) if pd.notna(p_val) else np.nan,
        })

    res_df = pd.DataFrame(rows)
    res_df = res_df.dropna(subset=["p_value"]).copy()

    # Benjamini-Hochberg False Discovery Rate (FDR) correction
    res_df["padj_bh"] = stats.false_discovery_control(res_df["p_value"], method="bh")

    # Sort most significant first
    res_df = res_df.sort_values(by=["padj_bh", "p_value"], ascending=[True, True]).reset_index(drop=True)
    return res_df


def main() -> None:
    print("Fetching joined data from PostgreSQL...")
    raw_df = fetch_joined_expression()
    print(f"Loaded {len(raw_df)} rows across {raw_df['ensg_id'].nunique()} genes.")

    print("Running differential expression statistics...")
    stats_df = run_analysis(raw_df)

    csv_path = OUTPUT_DIR / "mean_dge.csv"
    json_path = OUTPUT_DIR / "mean_dge.json"

    stats_df.to_csv(csv_path, index=False)
    stats_df.to_json(json_path, orient="records", indent=2)

    print(f"Saved {len(stats_df)} gene results to:")
    print(f"  - {csv_path}")
    print(f"  - {json_path}")
    print("\nTop 5 genes by FDR:")
    print(stats_df.head(5).to_string(index=False))


if __name__ == "__main__":
    main()