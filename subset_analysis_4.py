from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
SUBSET_OUTPUT = ROOT / "outputs" / "miraclib_baseline_samples.csv"
SUBSET_SUMMARY_OUTPUT = ROOT / "outputs" / "miraclib_baseline_summary.csv"

# Return samples matching all selected Part 4 filters.
SUBSET_QUERY = """
SELECT p.project_code project, sub.subject_code subject, sub.response,
       sub.sex, s.sample_code sample
FROM samples s
JOIN subjects sub ON s.subject_id = sub.id
JOIN projects p ON sub.project_id = p.id
WHERE sub.condition IN ({conditions})
  AND sub.treatment IN ({treatments})
  AND s.sample_type IN ({sample_types})
  AND s.time_from_treatment_start IN ({timepoints})
ORDER BY p.project_code, sub.subject_code;
"""


def as_list(value):
    return [value] if isinstance(value, (str, int)) else list(value)


def analyze_baseline_subset(
    conn,
    condition="melanoma",
    treatment="miraclib",
    sample_type="PBMC",
    timepoint=0,
    save_output=True,
):
    # Accept single values or dashboard multi-select values.
    filters = list(map(
        as_list, (condition, treatment, sample_type, timepoint)
    ))
    placeholders = [
        ", ".join(["?"] * len(values)) for values in filters
    ]
    query = SUBSET_QUERY.format(
        conditions=placeholders[0],
        treatments=placeholders[1],
        sample_types=placeholders[2],
        timepoints=placeholders[3],
    )
    params = [item for values in filters for item in values]
    df = pd.read_sql_query(query, conn, params=params)

    summaries = []
    for category, counted_column in (
        ("project", "sample"),
        ("response", "subject"),
        ("sex", "subject"),
    ):
        counts = df.groupby(category, as_index=False)[counted_column].nunique()
        counts.columns = ["value", "count"]
        counts.insert(0, "category", category)
        summaries.append(counts)
    summary_df = pd.concat(summaries, ignore_index=True)

    if save_output:
        df.to_csv(SUBSET_OUTPUT, index=False)
        summary_df.to_csv(SUBSET_SUMMARY_OUTPUT, index=False)
    return df, summary_df
