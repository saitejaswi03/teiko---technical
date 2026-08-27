#!/usr/bin/env python3
"""
Part 3: Statistical Analysis
=============================
Compares relative-frequency cell population profiles between miraclib
responders and non-responders, restricted to melanoma PBMC samples, in
order to find populations that might predict treatment response.

Outputs (written to output/):
    part3_responder_frequencies.csv   long-format table used for the test
    part3_stats_results.csv           per-population test statistics
    part3_boxplots.png                one boxplot per population

Statistical approach
---------------------
For each of the five populations we compare the responder vs. non-responder
relative-frequency distributions with:
  * Mann-Whitney U test (primary) - a non-parametric rank test that makes no
    normality assumption, appropriate given the modest per-group sample
    sizes and the fact that percentages are bounded and can be skewed.
  * Welch's t-test (secondary) - reported for comparison; does not assume
    equal variances.
Because five populations are tested simultaneously, raw p-values are also
adjusted for multiple testing with the Benjamini-Hochberg (FDR) procedure.
A population is called "significant" when its BH-adjusted Mann-Whitney
p-value < 0.05.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from analysis.db import get_connection
from analysis.part2_frequencies import compute_frequencies

OUTPUT_DIR = REPO_ROOT / "output"
FREQ_PATH = OUTPUT_DIR / "part3_responder_frequencies.csv"
STATS_PATH = OUTPUT_DIR / "part3_stats_results.csv"
PLOT_PATH = OUTPUT_DIR / "part3_boxplots.png"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


def get_responder_frequencies(conn) -> pd.DataFrame:
    """Melanoma, PBMC, miraclib samples with known response, with per-sample
    relative frequencies for every population (long format)."""
    freq = compute_frequencies(conn)

    meta = pd.read_sql_query(
        """
        SELECT s.sample_id AS sample, s.subject_id, s.sample_type,
               s.treatment, s.response, sub.condition
        FROM samples s
        JOIN subjects sub ON sub.subject_id = s.subject_id
        WHERE sub.condition = 'melanoma'
          AND s.sample_type = 'PBMC'
          AND s.treatment = 'miraclib'
          AND s.response IS NOT NULL
        """,
        conn,
    )

    merged = freq.merge(meta, on="sample", how="inner")
    return merged.sort_values(["population", "response", "sample"]).reset_index(drop=True)


def benjamini_hochberg(pvals: pd.Series) -> pd.Series:
    """Return BH-adjusted (FDR) q-values for a series of p-values."""
    n = len(pvals)
    order = pvals.sort_values().index
    ranked = pvals.loc[order]
    adj = ranked * n / pd.Series(range(1, n + 1), index=order)
    # enforce monotonicity (standard BH step-up correction)
    adj = adj[::-1].cummin()[::-1]
    return adj.clip(upper=1.0).reindex(pvals.index)


def run_tests(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for pop in POPULATIONS:
        sub = df[df.population == pop]
        yes = sub.loc[sub.response == "yes", "percentage"]
        no = sub.loc[sub.response == "no", "percentage"]

        u_stat, u_p = stats.mannwhitneyu(yes, no, alternative="two-sided")
        t_stat, t_p = stats.ttest_ind(yes, no, equal_var=False)

        rows.append(
            {
                "population": pop,
                "n_responders": len(yes),
                "n_non_responders": len(no),
                "mean_pct_responders": round(yes.mean(), 3),
                "mean_pct_non_responders": round(no.mean(), 3),
                "mannwhitney_u": round(u_stat, 3),
                "mannwhitney_p": u_p,
                "welch_t": round(t_stat, 3),
                "welch_p": t_p,
            }
        )

    result = pd.DataFrame(rows)
    result["mannwhitney_p_adj_bh"] = benjamini_hochberg(result["mannwhitney_p"])
    result["significant_bh_0.05"] = result["mannwhitney_p_adj_bh"] < 0.05
    return result.sort_values("mannwhitney_p_adj_bh").reset_index(drop=True)


def make_boxplots(df: pd.DataFrame, path: Path) -> None:
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(1, len(POPULATIONS), figsize=(4 * len(POPULATIONS), 5), sharey=False)

    for ax, pop in zip(axes, POPULATIONS):
        sub = df[df.population == pop]
        sns.boxplot(
            data=sub,
            x="response",
            y="percentage",
            order=["no", "yes"],
            hue="response",
            hue_order=["no", "yes"],
            palette={"no": "#d95f5f", "yes": "#4c8dae"},
            legend=False,
            ax=ax,
        )
        sns.stripplot(
            data=sub,
            x="response",
            y="percentage",
            order=["no", "yes"],
            color="black",
            alpha=0.35,
            size=3,
            ax=ax,
        )
        ax.set_title(pop)
        ax.set_xlabel("response")
        ax.set_ylabel("relative frequency (%)" if pop == POPULATIONS[0] else "")
        ax.set_xticks(range(2))
        ax.set_xticklabels(["non-responder", "responder"])

    fig.suptitle(
        "Melanoma PBMC samples on miraclib: relative frequency by response",
        y=1.03,
        fontsize=13,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    conn = get_connection()
    try:
        df = get_responder_frequencies(conn)
    finally:
        conn.close()

    if df.empty:
        raise RuntimeError(
            "No melanoma/PBMC/miraclib samples with a recorded response were found."
        )

    df.to_csv(FREQ_PATH, index=False)

    results = run_tests(df)
    results.to_csv(STATS_PATH, index=False)

    make_boxplots(df, PLOT_PATH)

    print(f"Wrote {FREQ_PATH.relative_to(REPO_ROOT)} ({len(df)} rows)")
    print(f"Wrote {STATS_PATH.relative_to(REPO_ROOT)}")
    print(f"Wrote {PLOT_PATH.relative_to(REPO_ROOT)}")
    print()
    print(results.to_string(index=False))


if __name__ == "__main__":
    main()
