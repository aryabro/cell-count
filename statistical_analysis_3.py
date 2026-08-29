from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.formula.api as smf

ROOT = Path(__file__).resolve().parent
STATS_OUTPUT = ROOT / "outputs" / "miraclib_response_stats.csv"
BOXPLOT_OUTPUT = ROOT / "outputs" / "miraclib_response_boxplot.png"

# Build the responder comparison from the selected sample groups.
RESPONSE_QUERY = """
SELECT s.sample_code sample, sub.subject_code subject, sub.response,
       cp.name population,
       100.0 * cc.count / SUM(cc.count) OVER (PARTITION BY s.id) percentage
FROM cell_counts cc
JOIN samples s ON cc.sample_id = s.id
JOIN subjects sub ON s.subject_id = sub.id
JOIN cell_populations cp ON cc.population_id = cp.id
WHERE sub.condition IN ({conditions})
  AND sub.treatment IN ({treatments})
  AND s.sample_type IN ({sample_types})
  {timepoint_filter}
  AND sub.response IN ('yes', 'no')
ORDER BY cp.id, s.sample_code;
"""


def as_list(value):
    return [value] if isinstance(value, (str, int)) else list(value)


def analyze_response_groups(
    conn,
    condition="melanoma",
    treatment="miraclib",
    sample_type="PBMC",
    timepoints=None,
    save_output=True,
):
    # Accept single values or dashboard multi-select values.
    conditions, treatments, sample_types = map(
        as_list, (condition, treatment, sample_type)
    )
    groups = [conditions, treatments, sample_types]
    params = [item for group in groups for item in group]
    timepoint_filter = ""

    if timepoints is not None:
        selected_timepoints = as_list(timepoints)
        timepoint_filter = (
            "AND s.time_from_treatment_start IN "
            f"({', '.join(['?'] * len(selected_timepoints))})"
        )
        params += selected_timepoints

    placeholders = [", ".join(["?"] * len(group)) for group in groups]
    query = RESPONSE_QUERY.format(
        conditions=placeholders[0],
        treatments=placeholders[1],
        sample_types=placeholders[2],
        timepoint_filter=timepoint_filter,
    )
    df = pd.read_sql_query(query, conn, params=params)
    df["response_bin"] = df["response"].map({"yes": 1, "no": 0})
    results = []

    for population, population_df in df.groupby("population", sort=False):
        responder = population_df["response"].eq("yes")
        non_responder = population_df["response"].eq("no")
        fitted = smf.mixedlm(
            "percentage ~ response_bin",
            population_df,
            groups=population_df["subject"],
        ).fit()
        effect = fitted.params["response_bin"]
        p_value = fitted.pvalues["response_bin"]
        results.append({
            "population": population,
            "responder_samples": responder.sum(),
            "non_responder_samples": non_responder.sum(),
            "responder_subjects": population_df.loc[
                responder, "subject"
            ].nunique(),
            "non_responder_subjects": population_df.loc[
                non_responder, "subject"
            ].nunique(),
            "response_effect": effect,
            "p_value": p_value,
            "significant": p_value < 0.05,
        })

    results_df = pd.DataFrame(results)
    if save_output:
        results_df.to_csv(STATS_OUTPUT, index=False)
    return df, results_df


def create_response_boxplot(
    df,
    condition="melanoma",
    treatment="miraclib",
    save_output=True,
):
    labels = [
        ", ".join(str(item).title() for item in as_list(value))
        for value in (condition, treatment)
    ]
    populations = df["population"].unique()
    figure, axes = plt.subplots(1, len(populations), figsize=(16, 5), sharey=True)

    for axis, population in zip(axes, populations):
        population_df = df[df["population"] == population]
        axis.boxplot(
            [
                population_df.loc[
                    population_df["response"].eq(response), "percentage"
                ]
                for response in ("yes", "no")
            ],
            tick_labels=["Responder", "Non-responder"],
        )
        axis.set_title(population)
        axis.tick_params(axis="x", rotation=45)

    axes[0].set_ylabel("Relative frequency (%)")
    figure.suptitle(
        f"Immune Cell Frequencies: {labels[0]} · {labels[1]} "
        "Responders vs Non-Responders"
    )
    plt.tight_layout()

    if save_output:
        figure.savefig(BOXPLOT_OUTPUT, dpi=150, bbox_inches="tight")
        plt.close(figure)
        return None
    return figure
