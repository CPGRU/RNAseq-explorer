import os
from pathlib import Path

import psycopg2
from dotenv import load_dotenv


# etl/src/etl/db_conn.py -> parents[3] is RNAseq-explorer/
ROOT_DIR = Path(__file__).resolve().parents[3]
ENV_PATH = ROOT_DIR / ".env"

load_dotenv(dotenv_path=ENV_PATH)


def get_pg_connection():
    missing = [
        var for var in ("DB_HOST", "DB_NAME", "DB_USER", "DB_PASSWORD")
        if not os.getenv(var)
    ]
    if missing:
        raise RuntimeError(
            f"Missing required environment variables in {ENV_PATH}: {', '.join(missing)}"
        )

    return psycopg2.connect(
        host=os.getenv("DB_HOST", "localhost"),
        port=os.getenv("DB_PORT", "5432"),
        dbname=os.getenv("DB_NAME"),
        user=os.getenv("DB_USER"),
        password=os.getenv("DB_PASSWORD"),
    )