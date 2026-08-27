#!/usr/bin/env python3
"""
Interactive dashboard for the Teiko technical assignment.

Run with:  make dashboard   (equivalent to `streamlit run dashboard/app.py`)

Reads directly from cell_counts.db (built by load_data.py) so it always
reflects the current database, and re-derives the Part 2/3/4 results live
rather than depending on files in output/.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from analysis.db import DB_PATH, get_connection
from analysis.part2_frequencies import compute_frequencies
from analysis.part3_stats import benjamini_hochberg
from analysis.part4_subset import (
    avg_b_cells_melanoma_male_responders_baseline,
    baseline_miraclib_melanoma_pbmc,
    samples_per_project,
    subjects_by_response,
    subjects_by_sex,
)

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

st.set_page_config(page_title="Teiko — Immune Cell Dashboard", layout="wide")


@st.cache_data(show_spinner=False)
def load_all():
    conn = get_connection()
    try:
        freq = compute_frequencies(conn)
        meta = pd.read_sql_query(
            """
            SELECT s.sample_id AS sample, s.subject_id, s.sample_type,
                   s.treatment, s.response, s.time_from_treatment_start,
                   sub.condition, sub.age, sub.sex, sub.project_id
            FROM samples s
            JOIN subjects sub ON sub.subject_id = s.subject_id
            """,
            conn,
        )
        cohort = baseline_miraclib_melanoma_pbmc(conn)
        avg_b = avg_b_cells_melanoma_male_responders_baseline(conn)
    finally:
        conn.close()
    full = freq.merge(meta, on="sample", how="left")
    return full, meta, cohort, avg_b


if not DB_PATH.exists():
    st.error(
        f"Database not found at `{DB_PATH.name}`. Run `python load_data.py` "
        "(or `make pipeline`) first."
    )
    st.stop()

full, meta, baseline_cohort, avg_b_cells = load_all()

st.title("Immune Cell Population Dashboard")
st.caption("Loblaw Bio — miraclib clinical trial · cell-count.csv")

tab2, tab3, tab4 = st.tabs(
    ["Part 2 · Frequencies", "Part 3 · Responder Analysis", "Part 4 · Subset Analysis"]
)

# ---------------------------------------------------------------- Part 2 ---
with tab2:
    st.subheader("Relative frequency of each cell population, per sample")
    st.write(
        "For every sample, total cell count is the sum across all five "
        "populations; percentage is each population's share of that total."
    )

    with st.expander("Filters", expanded=True):
        c1, c2, c3 = st.columns(3)
        projects = c1.multiselect("Project", sorted(meta.project_id.unique()))
        conditions = c2.multiselect("Condition", sorted(meta.condition.unique()))
        sample_types = c3.multiselect("Sample type", sorted(meta.sample_type.unique()))

    view = full.copy()
    if projects:
        view = view[view.project_id.isin(projects)]
    if conditions:
        view = view[view.condition.isin(conditions)]
    if sample_types:
        view = view[view.sample_type.isin(sample_types)]

    table_cols = ["sample", "total_count", "population", "count", "percentage"]
    st.dataframe(view[table_cols], use_container_width=True, height=350)
    st.download_button(
        "Download filtered table as CSV",
        view[table_cols].to_csv(index=False).encode(),
        file_name="cell_frequencies_filtered.csv",
        mime="text/csv",
    )

    st.markdown("**Population share distribution across filtered samples**")
    fig = px.box(view, x="population", y="percentage", points="outliers", color="population")
    fig.update_layout(showlegend=False, yaxis_title="relative frequency (%)")
    st.plotly_chart(fig, use_container_width=True)

# ---------------------------------------------------------------- Part 3 ---
with tab3:
    st.subheader("Miraclib responders vs. non-responders — melanoma, PBMC")
    sub = full[
        (full.condition == "melanoma")
        & (full.sample_type == "PBMC")
        & (full.treatment == "miraclib")
        & (full.response.notna())
    ]

    if sub.empty:
        st.warning("No matching samples in the database.")
    else:
        fig = px.box(
            sub,
            x="population",
            y="percentage",
            color="response",
            category_orders={"population": POPULATIONS, "response": ["no", "yes"]},
            color_discrete_map={"no": "#d95f5f", "yes": "#4c8dae"},
            points="all",
            labels={"percentage": "relative frequency (%)", "response": "response"},
        )
        fig.update_layout(boxmode="group")
        st.plotly_chart(fig, use_container_width=True)

        rows = []
        for pop in POPULATIONS:
            s = sub[sub.population == pop]
            yes = s.loc[s.response == "yes", "percentage"]
            no = s.loc[s.response == "no", "percentage"]
            if len(yes) < 2 or len(no) < 2:
                continue
            u, p = stats.mannwhitneyu(yes, no, alternative="two-sided")
            rows.append(
                {
                    "population": pop,
                    "n_responders": len(yes),
                    "n_non_responders": len(no),
                    "mean % (responders)": round(yes.mean(), 2),
                    "mean % (non-responders)": round(no.mean(), 2),
                    "Mann-Whitney p": p,
                }
            )
        results = pd.DataFrame(rows)
        results["BH-adjusted p"] = benjamini_hochberg(results["Mann-Whitney p"])
        results["significant (q<0.05)"] = results["BH-adjusted p"] < 0.05
        results = results.sort_values("BH-adjusted p").reset_index(drop=True)

        st.markdown("**Mann-Whitney U test, per population (Benjamini-Hochberg adjusted)**")
        st.dataframe(
            results.style.format(
                {
                    "Mann-Whitney p": "{:.4f}",
                    "BH-adjusted p": "{:.4f}",
                    "mean % (responders)": "{:.2f}",
                    "mean % (non-responders)": "{:.2f}",
                }
            ),
            use_container_width=True,
        )

        sig = results[results["significant (q<0.05)"]]
        if sig.empty:
            top = results.iloc[0]
            st.info(
                f"No population reaches significance after multiple-testing "
                f"correction (q < 0.05) in this dataset. `{top.population}` is "
                f"the most suggestive (BH-adjusted p = {top['BH-adjusted p']:.3f})."
            )
        else:
            st.success(
                "Significant populations (q < 0.05): "
                + ", ".join(sig.population.tolist())
            )

# ---------------------------------------------------------------- Part 4 ---
with tab4:
    st.subheader("Baseline cohort: melanoma, PBMC, miraclib, time = 0")

    c1, c2, c3 = st.columns(3)
    c1.metric("Samples", baseline_cohort["sample"].nunique())
    c2.metric("Subjects", baseline_cohort["subject_id"].nunique())
    c3.metric("Projects represented", baseline_cohort["project_id"].nunique())

    colA, colB, colC = st.columns(3)
    with colA:
        st.markdown("**Samples per project**")
        st.dataframe(samples_per_project(baseline_cohort), use_container_width=True)
    with colB:
        st.markdown("**Subjects by response**")
        st.dataframe(subjects_by_response(baseline_cohort), use_container_width=True)
    with colC:
        st.markdown("**Subjects by sex**")
        st.dataframe(subjects_by_sex(baseline_cohort), use_container_width=True)

    st.divider()
    st.markdown(
        "**Melanoma males, all sample & treatment types, responders, time = 0** "
        "— average B-cell count:"
    )
    st.metric("Average B cells", f"{avg_b_cells:.2f}")
