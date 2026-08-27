# Teiko Technical Assignment — Immune Cell Population Analysis

Analysis of `cell-count.csv` for Bob Loblaw's miraclib clinical trial: a
relational database, a reproducible analysis pipeline (Parts 1-4), and an
interactive dashboard.

## Quick start (GitHub Codespaces or any machine with Python 3.10+)

```bash
make setup       # installs dependencies from requirements.txt
make pipeline    # builds cell_counts.db and regenerates everything in output/
make dashboard   # starts the Streamlit dashboard (Codespaces will offer to
                  # forward the port and give you a browser link)
```

Running the three commands in order reproduces every table and figure in
this repository from `cell-count.csv` alone.

**Dashboard:** `make dashboard` runs `streamlit run dashboard/app.py` on
`localhost:8501`. In Codespaces, forward/open port 8501 to get a public
preview URL. *(If deployed to Streamlit Community Cloud, the live link is:
`<teiko---technical-2xmlk9xp7kemu7sx5f8u7s
.streamlit.app>` — see "Deploying the dashboard" below.)*

## Repository layout

```
load_data.py                 Part 1 — builds cell_counts.db from cell-count.csv
pipeline.py                  orchestrates load_data.py + all of analysis/ (what `make pipeline` runs)
analysis/
  db.py                      shared SQLite connection helper
  part2_frequencies.py       Part 2 — per-sample population frequency table
  part3_stats.py             Part 3 — responder vs non-responder stats + boxplot
  part4_subset.py            Part 4 — baseline cohort subset queries
dashboard/
  app.py                     Streamlit dashboard (Parts 2-4 in three tabs)
output/                      generated tables (.csv), figure (.png), summary (.txt)
cell-count.csv               input data
cell_counts.db               generated SQLite database (committed for convenience;
                              `make pipeline` regenerates it from scratch)
requirements.txt
Makefile
```

## Part 1 — Data management: schema design

`load_data.py` builds four normalized tables:

```
projects (project_id PK)
    |
subjects (subject_id PK, project_id FK, condition, age, sex)
    |
samples  (sample_id PK, subject_id FK, sample_type, treatment, response,
          time_from_treatment_start)
    |
cell_counts (sample_id FK, population, count)   -- one row per (sample, population)
             PK (sample_id, population)
```

**Why this shape, and not one wide table.** The source CSV is one row per
sample with the five cell populations as columns. That layout is convenient
to read but violates normal form twice over: subject-level attributes
(`condition`, `age`, `sex`) repeat across every sample from the same
subject, and the five population columns are really the same fact
("a count, for this sample, for this population") repeated five times as
different column names. Normalizing removes both redundancies:

- **`projects` / `subjects` / `samples`** separate three different grains
  of truth — a subject belongs to one project and contributes many samples
  over time; a sample belongs to one subject and carries the
  treatment/response/timepoint that describes *that draw*. Splitting them
  out means a subject's condition or age is stored once, not once per
  sample, so it can't drift out of sync and updates only touch one row.
- **`cell_counts` is long/tidy rather than wide.** Instead of five fixed
  columns, each (sample, population) pair is its own row. This is the
  design decision that matters most for scaling: adding a sixth immune
  population, or a panel with fifty populations, is an `INSERT`, not an
  `ALTER TABLE`. Every query in this project (Part 2's frequency table,
  Part 3's stats, Part 4's aggregates) is a `GROUP BY population`
  over this table, which is exactly the shape SQL aggregation is built for.

**How this scales to hundreds of projects, thousands of samples, and
new kinds of analysis.** The row counts change but the schema doesn't:

- Every foreign key here (`project_id`, `subject_id`, `sample_id`,
  `population`) is indexed (see the `CREATE INDEX` statements in
  `load_data.py`), so filtering — e.g. "melanoma, PBMC, miraclib,
  baseline" as in Part 4 — stays a fast indexed lookup rather than a full
  scan as the table grows into the millions of rows.
- New sample-level facts (a new assay readout, a QC flag, a batch ID) are
  new tables keyed on `sample_id`, not new columns bolted onto an
  ever-wider `samples` table. New subject-level facts (additional clinical
  metadata) similarly key on `subject_id`. This keeps each table's
  purpose narrow and avoids a monolithic table that every analysis has to
  wade through.
- SQLite is the right tool for this assignment's scale and the "just
  clone and run" grading requirement, but the schema itself is portable:
  the same DDL runs unchanged on Postgres, which is where you'd actually
  want to be once you have hundreds of projects, concurrent writers, or
  need row-level access control per project team.
- Because `cell_counts` is long-format, arbitrary "for population X,
  compare group A vs group B" analyses (Part 3's exact question) are one
  `GROUP BY population, response` away, regardless of how many
  populations exist.

## Part 2 — Data overview

`analysis/part2_frequencies.py` (invoked from `pipeline.py`, and
importable from the dashboard) computes, for every sample, the total cell
count across the five populations and each population's percentage of
that total. Output: `output/cell_frequencies.csv`, columns
`sample, total_count, population, count, percentage` (one row per
sample × population, 5 rows per sample).

## Part 3 — Statistical analysis

`analysis/part3_stats.py` restricts to **melanoma, PBMC, miraclib**
samples with a recorded response and compares each population's relative
frequency between responders and non-responders:

- **Mann-Whitney U test** (primary) — a non-parametric rank test, chosen
  because relative frequencies are bounded percentages that are not
  guaranteed to be normally distributed, and the test makes no such
  assumption.
- **Welch's t-test** (reported alongside, for readers who want the
  parametric comparison) — does not assume equal variances between groups.
- Because five populations are tested at once, raw p-values are corrected
  for multiple testing with the **Benjamini-Hochberg (FDR)** procedure;
  a population is called significant at **BH-adjusted p < 0.05**.

Outputs: `output/part3_boxplots.png` (one boxplot per population,
responders vs. non-responders, individual samples overlaid as points),
`output/part3_stats_results.csv` (test statistics per population), and
`output/part3_responder_frequencies.csv` (the underlying long-format data).

**Result on this dataset:** no population reaches significance after BH
correction; `cd4_t_cell` is the most suggestive (raw Mann-Whitney
p ≈ 0.013, BH-adjusted p ≈ 0.067) with responders showing a modestly
higher mean CD4 T-cell frequency. This is reported honestly rather than
rounded up to "significant" — see `output/part3_stats_results.csv` for the
full numbers, and the dashboard's Part 3 tab for the same result computed
live.

## Part 4 — Data subset analysis

`analysis/part4_subset.py` answers two independent questions:

**(A)** Melanoma, PBMC, miraclib, baseline (`time_from_treatment_start = 0`):
samples per project, subjects by response, subjects by sex.
See `output/part4_summary.txt` / `output/part4_samples_per_project.csv` /
`output/part4_responders_by_group.csv`.

**(B)** Melanoma males, *all* sample types and treatments, responders,
baseline (`time_from_treatment_start = 0`) — average B-cell count:

```
Average b_cell count: 10206.15
```

(computed directly from `cell_counts.db`; see
`avg_b_cells_melanoma_male_responders_baseline()` in
`analysis/part4_subset.py`, and `output/part4_summary.txt`.)

## Code structure rationale

- **`load_data.py` stays standalone and dependency-light** (stdlib
  `csv`/`sqlite3` only) so it satisfies the assignment's requirement to run
  as `python load_data.py` with zero setup beyond the CSV being present —
  it doesn't even need `pandas`.
- **`analysis/` is a plain importable package**, not a set of copy-pasted
  scripts: `part2_frequencies.compute_frequencies()` is called directly by
  `part3_stats.py` (Part 3 is Part 2's table, filtered and tested) and by
  the dashboard, so the frequency calculation is defined exactly once.
  Each `partN_*.py` module is also independently runnable
  (`python analysis/part3_stats.py`) for debugging one stage at a time.
- **`pipeline.py` is the single source of truth for "run everything"** —
  it's what `make pipeline` calls, and it's a thin sequence of calls into
  `load_data` and `analysis.*`, so there's exactly one place that defines
  pipeline order.
- **The dashboard queries the database directly** (via the same
  `analysis.*` functions used by the pipeline) rather than only reading the
  static CSVs in `output/`, so it always reflects the current
  `cell_counts.db` — re-run `make pipeline` with different/updated data and
  the dashboard picks it up on refresh.

## Deploying the dashboard (optional, for a shareable link)

`make dashboard` satisfies the grading requirement (a local server started
with one command). To also get a public URL to put in this README:

1. Push this repository to GitHub.
2. Go to [share.streamlit.io](https://share.streamlit.io), sign in, "New
   app", point it at this repo and `dashboard/app.py`.
3. Replace `<ADD_YOUR_DEPLOYED_URL_HERE>` above with the resulting URL.

## Requirements

See `requirements.txt`. Tested with Python 3.11/3.12.
