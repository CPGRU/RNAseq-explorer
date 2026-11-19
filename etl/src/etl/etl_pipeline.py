from pathlib import Path
import pandas as pd
from etl.db_conn import get_pg_connection

# Automatically find the project root (RNAseq-explorer) 
# Since etl_pipeline.py is in etl/src/etl/, we go up 3 levels to reach the root.
BASE_DIR = Path(__file__).resolve().parents[3]

# Build absolute paths to your CSV files
samples_csv_path = BASE_DIR / "data" / "archive" / "samples.csv"
expr_csv_path = BASE_DIR / "data" / "archive" / "gene_expression.csv"

# 1. Connect to your PostgreSQL database
conn = get_pg_connection()
cursor = conn.cursor()

# 2. Load and insert sample metadata CSV
# Ensure 'samples.csv' matches your actual filename
print("Importing samples metadata...")
samples_df = pd.read_csv(samples_csv_path)

for _, row in samples_df.iterrows():
    cursor.execute(
        """
        INSERT INTO samples (sample_id, dex, celltype, geo_id)
        VALUES (%s, %s, %s, %s)
        ON CONFLICT (sample_id) DO NOTHING;
        """,
        (row["sample_id"], row["dex"], row["celltype"], row["geo_id"])
    )
conn.commit()
print("Samples imported successfully!")

# 3. Load and insert gene expression CSV
# Ensure 'expression.csv' matches your actual filename
print("Importing gene expression data...")
expr_df = pd.read_csv(expr_csv_path)

# Melt the wide format into long format:
# - id_vars: columns to keep fixed (ensgene)
# - var_name: what to call the old column headers (sample_id)
# - value_name: what to call the cell values (expression_value)
melted_df = pd.melt(
    expr_df, 
    id_vars=["ensgene"], 
    var_name="sample_id", 
    value_name="expression_value"
)

# Iterate and insert the transformed long-format rows
for _, row in melted_df.iterrows():
    cursor.execute(
        """
        INSERT INTO gene_expression (sample_id, ensg_id, expression_value)
        VALUES (%s, %s, %s)
        ON CONFLICT (sample_id, ensg_id) DO NOTHING;
        """,
        (row["sample_id"], row["ensgene"], row["expression_value"])
    )

conn.commit()
print("Expressions imported successfully!")

# Close connection
cursor.close()
conn.close()
print("ETL process completed successfully!")