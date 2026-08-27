"""Small shared helper for connecting to cell_counts.db."""

from __future__ import annotations

import sqlite3
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DB_PATH = REPO_ROOT / "cell_counts.db"


def get_connection(db_path: Path | str = DB_PATH) -> sqlite3.Connection:
    db_path = Path(db_path)
    if not db_path.exists():
        raise FileNotFoundError(
            f"{db_path} not found. Run `python load_data.py` first to build "
            "the database from cell-count.csv."
        )
    return sqlite3.connect(db_path)
