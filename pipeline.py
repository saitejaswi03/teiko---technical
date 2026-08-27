#!/usr/bin/env python3
"""
pipeline.py
===========
Runs the entire data pipeline end to end, with no manual intervention:
  1. (Re)initialize the SQLite database and load cell-count.csv  (Part 1)
  2. Compute the per-sample population frequency table            (Part 2)
  3. Run the responder vs. non-responder statistical analysis
     and generate the boxplot figure                               (Part 3)
  4. Run the baseline-cohort subset queries                        (Part 4)

This is what `make pipeline` invokes. All outputs land in cell_counts.db
(repo root) and output/ (CSV tables, PNG figure, text summary).
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(REPO_ROOT))


def _banner(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def main() -> None:
    _banner("Part 1: load_data.py — initialize DB and load cell-count.csv")
    subprocess.run([sys.executable, str(REPO_ROOT / "load_data.py")], check=True)

    _banner("Part 2: cell population frequency table")
    from analysis import part2_frequencies

    part2_frequencies.main()

    _banner("Part 3: statistical analysis (responders vs non-responders)")
    from analysis import part3_stats

    part3_stats.main()

    _banner("Part 4: baseline cohort subset analysis")
    from analysis import part4_subset

    part4_subset.main()

    _banner("Pipeline complete — see output/ for all tables and figures")


if __name__ == "__main__":
    main()
