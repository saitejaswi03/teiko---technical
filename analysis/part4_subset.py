#!/usr/bin/env python3
"""
Part 4: Data Subset Analysis
==============================
Two independent queries against the database:

(A) Baseline miraclib melanoma PBMC cohort
    Filters to melanoma, PBMC, miraclib, time_from_treatment_start = 0, then
    reports:
      - number of samples per project
      - number of subjects who are responders / non-responders
      - number of subjects who are male / female

(B) Melanoma males, all sample types and treatments, at baseline (time=0)
    Average b_cell count among responders, to two decimal places.
    ("All sample and treatment types" means this query does NOT restrict to
    PBMC or to miraclib - only condition=melanoma, sex=male, response=yes,
    time_from_treatment_start=0.)

Outputs (written to output/):
    part4_samples_per_project.csv
    part4_responders_by_group.csv   (responder/non-responder and sex counts)
    part4_summary.txt               (human-readable summary incl. answer to B)
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from analysis.db import get_connection

OUTPUT_DIR = REPO_ROOT / "output"


def baseline_miraclib_melanoma_pbmc(conn) -> pd.DataFrame:
    """(A) melanoma, PBMC, miraclib, time_from_treatment_start = 0."""
    query = """
        SELECT s.sample_id AS sample, s.subject_id, s.response,
               sub.sex, sub.project_id
        FROM samples s
        JOIN subjects sub ON sub.subject_id = s.subject_id
        WHERE sub.condition = 'melanoma'
          AND s.sample_type = 'PBMC'
          AND s.treatment = 'miraclib'
          AND s.time_from_treatment_start = 0
    """
    return pd.read_sql_query(query, conn)


def samples_per_project(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby("project_id")["sample"]
        .nunique()
        .rename("n_samples")
        .reset_index()
        .sort_values("project_id")
    )


def subjects_by_response(df: pd.DataFrame) -> pd.DataFrame:
    subs = df.drop_duplicates("subject_id")
    return (
        subs.groupby("response")["subject_id"]
        .nunique()
        .rename("n_subjects")
        .reset_index()
    )


def subjects_by_sex(df: pd.DataFrame) -> pd.DataFrame:
    subs = df.drop_duplicates("subject_id")
    return (
        subs.groupby("sex")["subject_id"]
        .nunique()
        .rename("n_subjects")
        .reset_index()
    )


def avg_b_cells_melanoma_male_responders_baseline(conn) -> float:
    """(B) melanoma males, ALL sample types & treatments, responders, time=0.
    Average b_cell count, rounded to 2 decimals."""
    query = """
        SELECT cc.count
        FROM samples s
        JOIN subjects sub ON sub.subject_id = s.subject_id
        JOIN cell_counts cc ON cc.sample_id = s.sample_id
        WHERE sub.condition = 'melanoma'
          AND sub.sex = 'M'
          AND s.response = 'yes'
          AND s.time_from_treatment_start = 0
          AND cc.population = 'b_cell'
    """
    counts = pd.read_sql_query(query, conn)["count"]
    if counts.empty:
        raise RuntimeError("No matching samples found for the melanoma-male-responder query.")
    return round(float(counts.mean()), 2)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    conn = get_connection()
    try:
        cohort = baseline_miraclib_melanoma_pbmc(conn)
        by_project = samples_per_project(cohort)
        by_response = subjects_by_response(cohort)
        by_sex = subjects_by_sex(cohort)
        avg_b_cells = avg_b_cells_melanoma_male_responders_baseline(conn)
    finally:
        conn.close()

    by_project.to_csv(OUTPUT_DIR / "part4_samples_per_project.csv", index=False)

    responders_by_group = pd.concat(
        [
            by_response.rename(columns={"response": "group"}).assign(dimension="response"),
            by_sex.rename(columns={"sex": "group"}).assign(dimension="sex"),
        ],
        ignore_index=True,
    )[["dimension", "group", "n_subjects"]]
    responders_by_group.to_csv(OUTPUT_DIR / "part4_responders_by_group.csv", index=False)

    summary_lines = [
        "Part 4A: Melanoma, PBMC, miraclib, baseline (time_from_treatment_start = 0)",
        "-" * 78,
        f"Total samples in cohort: {cohort['sample'].nunique()}",
        f"Total subjects in cohort: {cohort['subject_id'].nunique()}",
        "",
        "Samples per project:",
        by_project.to_string(index=False),
        "",
        "Subjects by response:",
        by_response.to_string(index=False),
        "",
        "Subjects by sex:",
        by_sex.to_string(index=False),
        "",
        "Part 4B: Melanoma males, all sample & treatment types, responders, time=0",
        "-" * 78,
        f"Average b_cell count: {avg_b_cells:.2f}",
    ]
    summary = "\n".join(summary_lines)
    (OUTPUT_DIR / "part4_summary.txt").write_text(summary + "\n")

    print(summary)


if __name__ == "__main__":
    main()
