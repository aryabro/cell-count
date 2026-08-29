import sqlite3
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "data" / "cell-count.csv"
DB_PATH = ROOT / "cell-count.db"
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

# Keep repeated project, subject, and population data normalized.
SCHEMA = """
PRAGMA foreign_keys = ON;

CREATE TABLE projects (
    id INTEGER PRIMARY KEY,
    project_code TEXT NOT NULL UNIQUE
);

CREATE TABLE subjects (
    id INTEGER PRIMARY KEY,
    project_id INTEGER NOT NULL REFERENCES projects(id),
    subject_code TEXT NOT NULL,
    condition TEXT NOT NULL,
    age INTEGER CHECK (age IS NULL OR age >= 0),
    sex TEXT,
    treatment TEXT NOT NULL,
    response TEXT CHECK (response IS NULL OR response IN ('yes', 'no')),
    UNIQUE (project_id, subject_code)
);

CREATE TABLE samples (
    id INTEGER PRIMARY KEY,
    sample_code TEXT NOT NULL UNIQUE,
    subject_id INTEGER NOT NULL REFERENCES subjects(id),
    sample_type TEXT NOT NULL,
    time_from_treatment_start INTEGER NOT NULL
        CHECK (time_from_treatment_start >= 0)
);

CREATE TABLE cell_populations (
    id INTEGER PRIMARY KEY,
    name TEXT NOT NULL UNIQUE
);

CREATE TABLE cell_counts (
    sample_id INTEGER NOT NULL REFERENCES samples(id),
    population_id INTEGER NOT NULL REFERENCES cell_populations(id),
    count INTEGER NOT NULL CHECK (count >= 0),
    PRIMARY KEY (sample_id, population_id)
);
"""


def load_tables(conn, data):
    # Load parent tables before records that reference their IDs.
    projects = data[["project"]].drop_duplicates().rename(
        columns={"project": "project_code"}
    )
    projects.to_sql("projects", conn, if_exists="append", index=False)
    project_ids = dict(conn.execute("SELECT project_code, id FROM projects"))

    subject_columns = [
        "project", "subject", "condition", "age", "sex", "treatment", "response"
    ]
    subjects = data[subject_columns].drop_duplicates()
    subjects["project_id"] = subjects.pop("project").map(project_ids)
    subjects = subjects.rename(columns={"subject": "subject_code"})[
        [
            "project_id", "subject_code", "condition", "age",
            "sex", "treatment", "response",
        ]
    ]
    subjects.to_sql("subjects", conn, if_exists="append", index=False)
    subject_ids = {
        (project_id, code): subject_id
        for subject_id, project_id, code in conn.execute(
            "SELECT id, project_id, subject_code FROM subjects"
        )
    }

    project_id = data["project"].map(project_ids)
    samples = pd.DataFrame({
        "sample_code": data["sample"],
        "subject_id": [
            subject_ids[(pid, code)]
            for pid, code in zip(project_id, data["subject"])
        ],
        "sample_type": data["sample_type"],
        "time_from_treatment_start": data["time_from_treatment_start"],
    })
    samples.to_sql("samples", conn, if_exists="append", index=False)
    sample_ids = dict(conn.execute("SELECT sample_code, id FROM samples"))

    pd.DataFrame({"name": POPULATIONS}).to_sql(
        "cell_populations", conn, if_exists="append", index=False
    )
    population_ids = dict(
        conn.execute("SELECT name, id FROM cell_populations")
    )

    counts = data[["sample", *POPULATIONS]].melt(
        id_vars="sample", var_name="population", value_name="count"
    )
    counts["sample_id"] = counts.pop("sample").map(sample_ids)
    counts["population_id"] = counts.pop("population").map(population_ids)
    counts[["sample_id", "population_id", "count"]].to_sql(
        "cell_counts", conn, if_exists="append", index=False
    )


def main():
    data = pd.read_csv(CSV_PATH)

    # Rebuild from source to make every pipeline run reproducible.
    if DB_PATH.exists():
        DB_PATH.unlink()
    with sqlite3.connect(DB_PATH) as conn:
        conn.executescript(SCHEMA)
        load_tables(conn, data)


if __name__ == "__main__":
    main()
