import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FREQUENCY_OUTPUT = ROOT / "outputs" / "cell_frequencies.csv"


# Calculate one relative-frequency row per sample and population.
FREQUENCY_QUERY = """
SELECT
    s.sample_code AS sample,
    SUM(cc.count) OVER (
        PARTITION BY s.id
    ) AS total_count,
    cp.name AS population,
    cc.count AS count,
    100.0 * cc.count
        / SUM(cc.count) OVER (
            PARTITION BY s.id
        ) AS percentage
FROM cell_counts cc
JOIN samples s
    ON cc.sample_id = s.id
JOIN cell_populations cp
    ON cc.population_id = cp.id
ORDER BY s.sample_code, cp.id;
"""


def save_frequency_results(rows):
    # Preserve the column order requested in Part 2.
    with FREQUENCY_OUTPUT.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as file:
        writer = csv.writer(file)

        writer.writerow([
            "sample",
            "total_count",
            "population",
            "count",
            "percentage",
        ])

        for row in rows:
            writer.writerow([
                row["sample"],
                row["total_count"],
                row["population"],
                row["count"],
                round(row["percentage"], 4),
            ])
