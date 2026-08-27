#!/usr/bin/env python3
"""
Part 2: Data Overview
======================
Answers Bob's first question: "What is the frequency of each cell type in
each sample?"

Produces one row per (sample, population) with columns:
    sample, total_count, population, count, percentage

Can be run standalone (`python analysis/part2_frequencies.py`) or imported
and called via `compute_frequencies(conn)` from the pipeline / dashboard.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from analysis.db import get_connection
OUTPUT_PATH = REPO_ROOT / "output" / "cell_frequencies.csv"


def compute_frequencies(conn) -> pd.DataFrame:
    """Return the relative-frequency summary table described in Part 2."""
    counts = pd.read_sql_query(
        "SELECT sample_id AS sample, population, count FROM cell_counts", conn
    )
    totals = counts.groupby("sample")["count"].sum().rename("total_count")
    counts = counts.join(totals, on="sample")
    counts["percentage"] = (counts["count"] / counts["total_count"] * 100).round(4)

    counts = counts[["sample", "total_count", "population", "count", "percentage"]]
    return counts.sort_values(["sample", "population"]).reset_index(drop=True)


def main() -> None:
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    conn = get_connection()
    try:
        df = compute_frequencies(conn)
    finally:
        conn.close()

    df.to_csv(OUTPUT_PATH, index=False)
    print(f"Wrote {len(df)} rows to {OUTPUT_PATH.relative_to(REPO_ROOT)}")
    print(df.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
