import sqlite3
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

from statistical_analysis_3 import analyze_response_groups, create_response_boxplot
from subset_analysis_4 import analyze_baseline_subset

ROOT = Path(__file__).resolve().parent
DB = ROOT / "cell-count.db"
FREQUENCIES = ROOT / "outputs" / "cell_frequencies.csv"
FIELDS = ["condition", "treatment", "sample_type", "timepoint"]
LABELS = ["Conditions", "Treatments", "Sample types", "Timepoints"]


# Cache source data so dashboard reruns remain responsive.
@st.cache_data
def load_data():
    with sqlite3.connect(DB) as conn:
        counts = pd.read_sql_query("""
            SELECT 'Projects' label, COUNT(*) count FROM projects
            UNION ALL SELECT 'Subjects', COUNT(*) FROM subjects
            UNION ALL SELECT 'Samples', COUNT(*) FROM samples
            UNION ALL SELECT 'Populations', COUNT(*) FROM cell_populations
            UNION ALL SELECT 'Cell counts', COUNT(*) FROM cell_counts
        """, conn)
        response_options = pd.read_sql_query("""
            SELECT sub.condition, sub.treatment, s.sample_type,
                   s.time_from_treatment_start timepoint
            FROM subjects sub JOIN samples s ON s.subject_id = sub.id
            WHERE sub.response IN ('yes', 'no')
            GROUP BY sub.condition, sub.treatment, s.sample_type,
                     s.time_from_treatment_start
            HAVING COUNT(DISTINCT sub.response) = 2
        """, conn)
        subset_options = pd.read_sql_query("""
            SELECT DISTINCT sub.condition, sub.treatment, s.sample_type,
                   s.time_from_treatment_start timepoint
            FROM subjects sub JOIN samples s ON s.subject_id = sub.id
        """, conn)
    return pd.read_csv(FREQUENCIES), counts, response_options, subset_options


@st.cache_data(show_spinner=False)
def response_analysis(filters):
    with sqlite3.connect(DB) as conn:
        return analyze_response_groups(
            conn,
            condition=filters[0],
            treatment=filters[1],
            sample_type=filters[2],
            timepoints=filters[3],
            save_output=False,
        )


@st.cache_data(show_spinner=False)
def subset_analysis(filters):
    with sqlite3.connect(DB) as conn:
        return analyze_baseline_subset(
            conn,
            condition=filters[0],
            treatment=filters[1],
            sample_type=filters[2],
            timepoint=filters[3],
            save_output=False,
        )


def filter_controls(options, key, all_timepoints=False):
    # Each selection narrows the choices shown by the next filter.
    selected, remaining = [], options
    preferred = {
        "condition": "melanoma",
        "treatment": "miraclib",
        "sample_type": "PBMC",
        "timepoint": 0,
    }
    for column, field, label in zip(st.columns(4), FIELDS, LABELS):
        choices = sorted(remaining[field].dropna().unique().tolist())
        default = choices if field == "timepoint" and all_timepoints else [
            preferred[field] if preferred[field] in choices else choices[0]
        ] if choices else []
        values = column.multiselect(
            label, choices, default=default, key=f"{key}_{field}"
        )
        selected.append(tuple(values))
        remaining = remaining[remaining[field].isin(values)]
    return tuple(selected)


def pie_chart(data, title):
    st.subheader(title)
    if data.empty:
        st.info("No data available.")
        return
    total = data["count"].sum()
    figure, axis = plt.subplots(figsize=(5, 4))
    axis.pie(
        data["count"],
        labels=data["value"],
        autopct=lambda percentage: (
            f"{round(percentage * total / 100):,} "
            f"({percentage:.1f}%)"
        ),
        startangle=90,
    )
    axis.axis("equal")
    st.pyplot(figure)
    plt.close(figure)


st.set_page_config(
    page_title="Loblaw Bio Clinical Trial",
    page_icon="🔬",
    layout="wide",
)
st.title("Loblaw Bio Clinical Trial")

if not DB.exists() or not FREQUENCIES.exists():
    st.error("Run `make pipeline` before starting the dashboard.")
    st.stop()

frequency_df, counts_df, response_options, subset_options = load_data()

# Organize the dashboard around the four assignment sections.
part_1, part_2, part_3, part_4 = st.tabs([
    "Part 1 · Data Management",
    "Part 2 · Initial Analysis",
    "Part 3 · Statistical Analysis",
    "Part 4 · Subset Analysis",
])

with part_1:
    st.header("Relational database overview")
    columns = st.columns(len(counts_df))
    for column, row in zip(columns, counts_df.itertuples()):
        column.metric(row.label, f"{row.count:,}")
    st.code(
        "projects 1 ── many subjects 1 ── many samples\n"
        "samples 1 ── many cell_counts many ── 1 cell_populations"
    )
    st.markdown("""
    - **projects:** clinical projects
    - **subjects:** patient details, treatment, and response
    - **samples:** sample type and treatment timepoint
    - **cell_populations:** immune population names
    - **cell_counts:** population counts for each sample
    """)

with part_2:
    st.header("Cell population frequencies")
    sample = st.selectbox("Choose a sample", frequency_df["sample"].unique())
    sample_df = frequency_df[frequency_df["sample"] == sample]
    st.metric("Total cells", f"{int(sample_df['total_count'].iloc[0]):,}")
    chart, table = st.columns([2, 1])
    with chart:
        figure, axis = plt.subplots(figsize=(7, 5))
        axis.pie(
            sample_df["percentage"],
            labels=sample_df["population"],
            autopct="%1.1f%%",
            startangle=90,
        )
        st.pyplot(figure)
        plt.close(figure)
    table.dataframe(sample_df, hide_index=True, width="stretch")
    with st.expander("Complete frequency table"):
        st.dataframe(frequency_df, hide_index=True, width="stretch")
        st.download_button(
            "Download CSV",
            frequency_df.to_csv(index=False),
            "cell_frequencies.csv",
        )

with part_3:
    st.header("Responders versus non-responders")
    st.caption(
        "Results are pooled across selected groups. "
        "The model accounts for repeated samples from each subject."
    )
    filters = filter_controls(response_options, "response", all_timepoints=True)
    if not all(filters):
        st.warning("Select at least one value in every filter.")
    else:
        with st.spinner("Calculating statistics and boxplots..."):
            response_df, stats_df = response_analysis(*filters)
            figure = create_response_boxplot(
                response_df,
                condition=filters[0],
                treatment=filters[1],
                save_output=False,
            )
        significant = stats_df.loc[stats_df["significant"], "population"]
        message = (
            f"Significant difference: {', '.join(significant)}."
            if len(significant)
            else "No significant differences at p < 0.05."
        )
        st.info(message)
        st.dataframe(stats_df, hide_index=True, width="stretch")
        st.pyplot(figure)
        plt.close(figure)

with part_4:
    st.header("Filtered sample subset")
    filters = filter_controls(subset_options, "subset")
    if not all(filters):
        st.warning("Select at least one value in every filter.")
    else:
        with st.spinner("Loading sample subset..."):
            subset_df, summary_df = subset_analysis(*filters)
        samples, subjects = st.columns(2)
        samples.metric("Matching samples", f"{subset_df['sample'].nunique():,}")
        subjects.metric("Subjects", f"{subset_df['subject'].nunique():,}")
        charts = st.columns(3)
        sections = [
            ("project", "Samples per project"),
            ("response", "Subjects by response"),
            ("sex", "Subjects by sex"),
        ]
        for column, (category, title) in zip(charts, sections):
            with column:
                pie_chart(summary_df[summary_df["category"] == category], title)
        st.dataframe(subset_df, hide_index=True, width="stretch")
        st.download_button(
            "Download CSV",
            subset_df.to_csv(index=False),
            "filtered_samples.csv",
        )
