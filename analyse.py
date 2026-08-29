import sqlite3
from pathlib import Path

from initial_analysis_2 import (
    FREQUENCY_QUERY,
    save_frequency_results,
)
from statistical_analysis_3 import (
    analyze_response_groups,
    create_response_boxplot,
)
from subset_analysis_4 import (
    analyze_baseline_subset,
)

ROOT = Path(__file__).resolve().parent

# ---------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------

def main():
    (ROOT / "outputs").mkdir(exist_ok=True)

    # Run Parts 2–4 against the database created by load_data.py.
    with sqlite3.connect(ROOT / "cell-count.db") as conn:
        conn.row_factory = sqlite3.Row

        # Part 2
        frequency_rows = conn.execute(FREQUENCY_QUERY).fetchall()
        save_frequency_results(frequency_rows)

        # Part 3
        response_df, _ = analyze_response_groups(conn)

        # Part 4
        analyze_baseline_subset(conn)

    create_response_boxplot(response_df)


if __name__ == "__main__":
    main()