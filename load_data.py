#!/usr/bin/env python3
"""
load_data.py
============
Part 1 of the Teiko technical assignment: Data Management.

Initializes a SQLite database (`cell_counts.db`, created in the repository
root) with a normalized relational schema, then loads every row of
`cell-count.csv` into it.

Usage
-----
    python load_data.py

No command-line arguments are required. The script is safe to re-run: it
drops and recreates the schema each time so the database always reflects
the current contents of cell-count.csv.

Schema
------
projects        (project_id)
subjects        (subject_id, project_id, condition, age, sex)
samples         (sample_id, subject_id, sample_type, treatment, response,
                 time_from_treatment_start)
cell_counts     (sample_id, population, count)   -- long/tidy format,
                 one row per (sample, immune population) pair.

See README.md for the full rationale behind this design.
"""

from __future__ import annotations

import csv
import sqlite3
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
CSV_PATH = REPO_ROOT / "cell-count.csv"
DB_PATH = REPO_ROOT / "cell_counts.db"

# The five immune cell population columns present in cell-count.csv.
# Keeping this as an explicit list (rather than inferring it) makes the
# long/tidy transform below unambiguous and easy to extend.
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

SCHEMA = """
DROP TABLE IF EXISTS cell_counts;
DROP TABLE IF EXISTS samples;
DROP TABLE IF EXISTS subjects;
DROP TABLE IF EXISTS projects;

CREATE TABLE projects (
    project_id      TEXT PRIMARY KEY
);

CREATE TABLE subjects (
    subject_id      TEXT PRIMARY KEY,
    project_id      TEXT NOT NULL REFERENCES projects(project_id),
    condition       TEXT NOT NULL,      -- e.g. melanoma, carcinoma, healthy
    age             INTEGER,
    sex             TEXT
);

CREATE TABLE samples (
    sample_id                  TEXT PRIMARY KEY,
    subject_id                 TEXT NOT NULL REFERENCES subjects(subject_id),
    sample_type                TEXT NOT NULL,   -- PBMC, WB, ...
    treatment                  TEXT NOT NULL,   -- miraclib, phauximab, none
    response                   TEXT,            -- 'yes' / 'no' / NULL
    time_from_treatment_start  INTEGER
);

CREATE TABLE cell_counts (
    sample_id       TEXT NOT NULL REFERENCES samples(sample_id),
    population      TEXT NOT NULL,     -- b_cell, cd8_t_cell, cd4_t_cell, nk_cell, monocyte
    count           INTEGER NOT NULL,
    PRIMARY KEY (sample_id, population)
);

CREATE INDEX idx_subjects_project   ON subjects(project_id);
CREATE INDEX idx_samples_subject    ON samples(subject_id);
CREATE INDEX idx_samples_filters    ON samples(sample_type, treatment, response, time_from_treatment_start);
CREATE INDEX idx_cellcounts_sample  ON cell_counts(sample_id);
CREATE INDEX idx_cellcounts_pop     ON cell_counts(population);
"""


def init_db(conn: sqlite3.Connection) -> None:
    """Create (or recreate) the schema described in SCHEMA."""
    conn.executescript(SCHEMA)


def load_csv(conn: sqlite3.Connection, csv_path: Path) -> int:
    """Read cell-count.csv and populate all four tables. Returns row count."""
    if not csv_path.exists():
        raise FileNotFoundError(
            f"Could not find {csv_path}. Make sure cell-count.csv is in the "
            "repository root before running load_data.py."
        )

    projects_seen: set[str] = set()
    subjects_seen: set[str] = set()
    n_rows = 0

    cur = conn.cursor()
    with csv_path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            project_id = row["project"]
            subject_id = row["subject"]
            sample_id = row["sample"]

            if project_id not in projects_seen:
                cur.execute(
                    "INSERT OR IGNORE INTO projects (project_id) VALUES (?)",
                    (project_id,),
                )
                projects_seen.add(project_id)

            if subject_id not in subjects_seen:
                cur.execute(
                    """INSERT OR IGNORE INTO subjects
                       (subject_id, project_id, condition, age, sex)
                       VALUES (?, ?, ?, ?, ?)""",
                    (subject_id, project_id, row["condition"], int(row["age"]), row["sex"]),
                )
                subjects_seen.add(subject_id)

            response = row["response"] if row["response"] not in ("", None) else None
            cur.execute(
                """INSERT INTO samples
                   (sample_id, subject_id, sample_type, treatment, response,
                    time_from_treatment_start)
                   VALUES (?, ?, ?, ?, ?, ?)""",
                (
                    sample_id,
                    subject_id,
                    row["sample_type"],
                    row["treatment"],
                    response,
                    int(row["time_from_treatment_start"]),
                ),
            )

            for pop in POPULATIONS:
                cur.execute(
                    """INSERT INTO cell_counts (sample_id, population, count)
                       VALUES (?, ?, ?)""",
                    (sample_id, pop, int(row[pop])),
                )

            n_rows += 1

    conn.commit()
    return n_rows


def main() -> None:
    # Start from a clean database file every run so `python load_data.py`
    # is idempotent and always matches the current CSV.
    if DB_PATH.exists():
        DB_PATH.unlink()

    conn = sqlite3.connect(DB_PATH)
    try:
        init_db(conn)
        n_rows = load_csv(conn, CSV_PATH)
        print(f"Loaded {n_rows} rows from {CSV_PATH.name} into {DB_PATH.name}")

        n_projects = conn.execute("SELECT COUNT(*) FROM projects").fetchone()[0]
        n_subjects = conn.execute("SELECT COUNT(*) FROM subjects").fetchone()[0]
        n_samples = conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
        n_counts = conn.execute("SELECT COUNT(*) FROM cell_counts").fetchone()[0]
        print(
            f"  projects={n_projects}  subjects={n_subjects}  "
            f"samples={n_samples}  cell_counts={n_counts}"
        )
    finally:
        conn.close()


if __name__ == "__main__":
    main()
