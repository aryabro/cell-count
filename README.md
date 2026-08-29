# Immune Cell Trial Analysis

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/aryabro/cell-count)

## Database schema

The source CSV contains project, subject, sample, treatment, response, and immune cell count information in a single wide table. I normalized this data into five SQLite tables:

```text
projects
   |
   | 1:N
   v
subjects
   |
   | 1:N
   v
samples
   |
   | 1:N
   v
cell_counts
   |
   | N:1
   v
cell_populations
```

### `projects`

Stores each clinical project once.

```text
id
project_code
```

`project_code` is unique, while the integer `id` is used internally as the primary key.

A separate project table avoids repeating project identifiers throughout the database and gives the model a natural place for additional project-level metadata in the future.

### `subjects`

Stores attributes that are constant for a subject in the provided dataset.

```text
id
project_id
subject_code
condition
age
sex
treatment
response
```

Subjects belong to projects through `project_id`.

The combination of:

```text
(project_id, subject_code)
```

is unique rather than assuming that a subject identifier must be globally unique across every project.

Treatment and response are stored at the subject level because they remain constant across all samples for a subject in the supplied dataset. If future datasets contained treatment changes or longitudinal response assessments, these could instead be represented by separate treatment or response-event tables.

### `samples`

Stores sample-specific information.

```text
id
sample_code
subject_id
sample_type
time_from_treatment_start
```

Each sample belongs to one subject.

This separates subject metadata such as age and treatment from sample metadata such as sample type and collection timepoint.

### `cell_populations`

Stores immune cell population names.

```text
id
name
```

The provided dataset contains five populations:

```text
b_cell
cd8_t_cell
cd4_t_cell
nk_cell
monocyte
```

The database is not limited to these five populations.

### `cell_counts`

Stores the measured count for a population in a sample.

```text
sample_id
population_id
count
```

The primary key is:

```text
(sample_id, population_id)
```

so a sample can only have one stored count for a given population.

### Why cell counts are stored in long format?

The input CSV stores populations as separate columns:

```text
sample | b_cell | cd8_t_cell | cd4_t_cell | nk_cell | monocyte
```

During ingestion, these columns are converted into rows:

```text
sample_id | population_id | count
1         | 1             | 10908
1         | 2             | 24440
1         | 3             | 20491
...
```

This was the main extensibility decision in the schema.

If another immune population is introduced, it can be added to `cell_populations` and stored in `cell_counts` without adding a new database column or changing the analytical schema.

This becomes increasingly useful as the number of measured populations grows.

### Scaling the schema

The design is intended to work without structural changes if the dataset grows to hundreds of projects and thousands of samples.

Project and subject metadata is stored once rather than repeated for every sample. Likewise, population names are stored once instead of being repeated as free-form text in every cell-count record.

The large analytical table is therefore `cell_counts`, whose rows contain compact foreign keys and the measured count.

The schema also keeps raw measurements separate from derived analytical values. For example, relative frequencies are calculated from cell counts rather than permanently stored in the database, preventing derived values from becoming inconsistent with their source data.

Additional indexes could be introduced later based on real query patterns if the dataset grows substantially, without changing the relational model.

## Analysis

### Part 2 - Cell population frequencies

For each sample, the total number of cells is calculated by summing the counts of all measured populations.

The relative frequency of each population is then:

```text
percentage = population count / total sample count * 100
```

The resulting table contains:

```text
sample
total_count
population
count
percentage
```

The calculation is performed directly from the normalized `cell_counts` table using a SQL window function.

The output is written to:

```text
outputs/cell_frequencies.csv
```

This produces one row for each sample/population combination.

## Statistical analysis

Part 3 compares immune cell relative frequencies between responders and non-responders among:

```text
condition = melanoma
treatment = miraclib
sample_type = PBMC
response = yes or no
```

Subjects in the dataset *can contribute samples from multiple treatment timepoints*. These measurements are therefore not fully independent observations.

For this reason, I use a linear mixed-effects model for each immune population:

```text
percentage ~ response
```

with subject as a random effect.

This allows all qualifying samples to contribute to the analysis while accounting for repeated measurements from the same subject.

Response is encoded as:

```text
no  = 0
yes = 1
```

The response coefficient therefore represents the estimated percentage-point difference between responders and non-responders.

For each population, the analysis reports:

```text
responder sample count
non-responder sample count
responder subject count
non-responder subject count
response effect (percentage points)
p-value
significant
```

A significance threshold of:

```text
p < 0.05
```

is used.

The statistical results are saved to:

```text
outputs/miraclib_response_stats.csv
```

and the responder/non-responder distributions are visualized as boxplots in:

```text
outputs/miraclib_response_boxplot.png
```

## Part 4 - Baseline subset analysis

Part 4 selects samples matching:

```text
condition = melanoma
treatment = miraclib
sample_type = PBMC
time_from_treatment_start = 0
```

For the resulting subset, the analysis reports:

- number of samples from each project
- number of unique subjects by response
- number of unique subjects by sex

The qualifying samples are written to:

```text
outputs/miraclib_baseline_samples.csv
```

and the grouped summary is written to:

```text
outputs/miraclib_baseline_summary.csv
```

## Code structure

The code is separated by responsibility so that data ingestion, individual analyses, orchestration, and presentation remain independent.

### `load_data.py`

Creates the SQLite database and loads `Data/cell-count.csv`.

The main `load_tables()` function:

1. loads unique projects
2. loads unique subjects
3. creates samples and connects them to subjects
4. loads immune population definitions
5. converts the wide population columns into long-format cell-count records

The database is rebuilt from the source CSV each time the script runs so repeated pipeline executions produce the same result.

The generated database is:

```text
cell-count.db
```

in the repository root.

### `initial_analysis_2.py`

Contains the SQL query for Part 2 and the function that writes the resulting frequency summary.

`FREQUENCY_QUERY` calculates sample totals and relative frequencies directly from the database.

`save_frequency_results()` writes the required output columns to `cell_frequencies.csv`.

### `statistical_analysis_3.py`

Contains the responder versus non-responder analysis.

`analyze_response_groups()`:

- queries the required clinical subset
- calculates relative frequencies
- fits one mixed-effects model per immune population
- reports response effect (percentage points) sizes and p-values

`create_response_boxplot()` creates the responder/non-responder population boxplots.

The functions also accept filter values so the same analysis logic can be reused by the interactive dashboard.

### `subset_analysis_4.py`

Contains the Part 4 filtering and summary logic.

`analyze_baseline_subset()` queries the requested sample subset and calculates:

- samples by project
- subjects by response
- subjects by sex

The function also accepts alternative filter selections so the dashboard can explore related subsets without duplicating the query logic.

### `analyse.py`

Acts as the analysis pipeline entry point.

It runs Parts 2–4 sequentially against the database created by `load_data.py`.

Keeping orchestration here allows the individual analysis modules to remain independently reusable by both the batch pipeline and the dashboard.

### `dashboard.py`

Contains the Streamlit interface.

The dashboard is organized around the four parts of the assignment:

```text
Part 1 · Data Management
Part 2 · Initial Analysis
Part 3 · Statistical Analysis
Part 4 · Subset Analysis
```

It reuses the same analysis functions used by the pipeline rather than implementing separate statistical logic.

The dashboard also provides filters for exploring related treatment groups, sample types, conditions, and timepoints.

## Running the project

The project is designed to run from the repository root in GitHub Codespaces.

### 1. Install dependencies

```bash
make setup
```

This installs the dependencies listed in `requirements.txt`.

### 2. Run the complete pipeline

```bash
make pipeline
```

This runs:

```text
load_data.py
    ↓
cell-count.db
    ↓
analyse.py
    ↓
Parts 2–4 outputs
```

The generated outputs are stored under:

```text
outputs/
```

### 3. Start the dashboard

```bash
make dashboard
```

This starts the Streamlit application.

If using GitHub Codespaces, open the forwarded Streamlit port when prompted.

## Dashboard

Interactive dashboard:

[Open the dashboard](TODO)

## References

The implementation and analysis choices were informed by the following documentation and background reading:

- [SQLite Foreign Key Support](https://www.sqlite.org/foreignkeys.html) - reference for relational integrity and SQLite foreign-key behavior.
- [SQLite Window Functions](https://www.sqlite.org/windowfunctions.html) - used for the per-sample total and relative-frequency calculations in Part 2.
- [pandas `melt`](https://pandas.pydata.org/docs/reference/api/pandas.melt.html) - used in `load_data.py` to reshape the five cell-population columns from wide to long format before loading `cell_counts`.
- [statsmodels `MixedLM`](https://www.statsmodels.org/stable/generated/statsmodels.regression.mixed_linear_model.MixedLM.html) - API reference for the linear mixed-effects model used in Part 3.
- [Repeated Measures Designs and Analysis of Longitudinal Data](https://pmc.ncbi.nlm.nih.gov/articles/PMC6072386/) - background on why repeated observations from the same subject are correlated and should not be treated as independent measurements.
- [Streamlit caching documentation](https://docs.streamlit.io/develop/concepts/architecture/caching) - reference for caching database-derived data and analysis results in the interactive dashboard.
- [Tidy data paper](https://www.jstatsoft.org/article/view/v059i10) - Paper advising structuring data so each observation is a row and each variable is a column.
