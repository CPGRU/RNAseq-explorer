import pandas as pd
from etl.db_conn import get_pg_connection

# Connect to PostgreSQL
conn = get_pg_connection()

# Load joined data directly into a pandas DataFrame
query = """
    SELECT s.sample_id, s.dex, s.celltype, g.ensg_id, g.expression_value
    FROM samples s
    JOIN gene_expression g ON s.sample_id = g.sample_id;
"""
df = pd.read_sql(query, conn)
conn.close()

# Now you can use pandas/scipy/numpy for grouping, filtering, or statistical analysis
print( df.head(10))
print(df.groupby('dex')['expression_value'].mean())
print("test")