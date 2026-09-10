# Improving Parent Identification in Plant Breeding

*This repository accompanies the manuscript currently under review.*

**Improving parent identification in plant breeding using combined similarity and likelihood metrics**

Authors: *Alexandra Nirsha¹,  Martin Spanoghe², John Riviere², Deborah Lanterbecq¹,²,³*

¹ Centre pour l’Agronomie et l’Agro-Industrie de la Province de Hainaut (asbl CARAH), 7800 Ath, Belgium

² Laboratoire de Biotechnologie et Biologie Appliquée, Haute Ecole Provinciale de Hainaut-CONDORCET, 7000 Mons, Belgium

³ Hainaut Analyses (HA), 7000 Mons, Belgium

Journal: *Discover Plants* (under review)

---

# Overview

This repository provides a complete workflow for evaluating parent–offspring relationships from binary molecular marker data.

The implemented pipeline combines:

* Jaccard similarity analysis;
* LOD score computation;
* comparative analysis of both approaches;
* graphical visualization of the results.

The workflow reproduces the analyses presented in the accompanying manuscript.

## Contents

* anonymized example dataset.
* pedigree preparation functions;
* Jaccard similarity analysis;
* LOD score analysis;
* combined Jaccard–LOD comparison;
* plotting utilities;
* scripts reproducing the manuscript results;

## Workflow

The analysis follows the workflow below:

1. Load the project and required dependencies.
2. Load and prepare the binary data and pedigree information.
3. Define the relevant pedigree subsets.
4. Perform Jaccard similarity analysis.
5. Perform LOD score analysis.
6. Combine the results from both approaches.
7. Generate figures reported in the manuscript.

```text
Input binary data
        +
Pedigree metadata
        │
        ▼
Data & pedigree preparation
        │
        ▼
Similarity and likelihood computation
        │
      ┌─┴─┐
      ▼   ▼
  Jaccard  LOD
  analysis analysis
      │   │
      └─┬─┘
        ▼
 Combined analysis
        │
        ▼
Tables and figures
```

---

# Installation

## Clone the repository

```bash
git clone https://github.com/ale-nir/parent-identification-plant-breeding.git
```

Alternatively, download the repository as a ZIP archive from GitHub and extract it.

---

## Repository structure

```text
.
├── README.md
├── LICENSE
├── setup.R
│
├── data/
│	├── bin.csv
│   └── data.csv
│
├── docs/
│   ├── 01-pedigree_preparation.md
│   ├── 02-jaccard_analysis.md
│   ├── 03-lod_analysis.md
│	├── 04-combined_analysis.md
│	├── 05-plotting.md
│   └── 06_tables.md
│
├── R/
│	├── data_io.R
│	├── pedigree_processing.R
│	├── jaccard_analysis.R
│	├── lod_analysis.R
│	├── combined_analysis.R
│	├── parentage_utils.R
│	└── plotting_utils.R
│
├── src/
│	├── jaccard_cpp.cpp
│   └── lod_cpp.cpp
│
└── examples/
    ├── reproduce_analysis_workflow.R
    ├── reproduce_manuscript_figures.R
    └── reproduce_manuscript_tables.R
```

---

## Open the project

Open the project in RStudio (recommended), or set the working directory to the repository root.

```r
setwd("path/to/repository")
```

---

## Quick start

Source the setup script:

```r
source("setup.R")
```

This script automatically:

* installs missing packages, if necessary;
* loads all required packages and functions;
* compiles the C++ source files;
* loads the example dataset.

---

# Minimal reproducible workflow

For a complete reproduction of the analysis without running each step manually, use:

```r
source("examples/reproduce_analysis_workflow.R")
```

This script runs the complete analytical workflow, from pedigree preparation through the combined LOD–Jaccard analysis, and generates all main objects required for the downstream figures and tables.

The workflow generates:

* pedigree subsets: `varieties_both_parents` and `varieties_any_parent`;
* Jaccard results: `J`, `jaccard_both`, `jaccard_any`, and the corresponding structural-analysis results;
* LOD results: `L`, `lod_both`, and `lod_any`;
* combined LOD–Jaccard results: `combined_both`, `combined_any`, `combined_both_iqr`, and `combined_any_iqr`.

The detailed workflow described below is optional and can be followed to reproduce or inspect each analysis step individually. If the complete workflow has already been executed using `examples/reproduce_analysis_workflow.R`, you can proceed directly to Steps 5 and 6 to generate the manuscript figures and tables.

# Reproduce the Analysis

The following sections describe the analysis workflow step by step. Each step provides instructions and identifies the main objects generated for use in downstream analyses.

## Step 1 – Prepare the data

The main datasets are loaded when `setup.R` is sourced. The data must then be processed according to the instructions below:

[Pedigree preparation](docs/01-pedigree_preparation.md)

After completing this step, two main objects should be available, corresponding to the biparental and single-parent subsets:

```r
varieties_both_parents
varieties_any_parent
```

These objects should be retained for downstream analyses.

---

## Step 2 – Compute Jaccard similarity for both subsets

Follow the instructions below to compute Jaccard similarity scores for the two pedigree subsets:

[Jaccard analysis](docs/02-jaccard_analysis.md)

This step first computes a complete pairwise Jaccard similarity matrix and then uses it for the parent-retrieval and structural analyses.

The following main objects are generated:

### 1. Pairwise Jaccard similarity matrix

```r
J
```

`J` contains the pairwise Jaccard similarity scores for all varieties in the dataset. Each value represents the Jaccard similarity between a pair of varieties and is used as the input for all subsequent Jaccard analyses.

The matrix is computed once and reused throughout the workflow.

### 2. Jaccard parent-retrieval results

```r
jaccard_both
jaccard_any
```

These objects contain the parent-retrieval results for different top-*k* values (*k* = 1, 5, 10, and 20).

To access the results for a specific top-*k* value, use, for example:

```r
jaccard_both$Top10
```

The returned object contains:

* detailed parent-retrieval results;
* recovery statistics;
* summary tables.

### 3. Relationship-structure results

The following objects contain information on the relationship structure used to evaluate the Jaccard similarity distributions:

```r
jaccard_struct_both
jaccard_struct_any
```

The relationship groups are:

* **Parent** – documented parent–offspring relationships;
* **Sibling** – varieties sharing at least one documented parent;
* **Random** – randomly selected unrelated varieties.

### 4. Similarity-distribution results

The following objects contain detailed Jaccard similarity distributions for the different relationship groups:

```r
results_both
results_any
```

Each object contains:

* `summary` – descriptive statistics for each relationship group;
* `kruskal` – Kruskal–Wallis test results;
* `dunn` – Dunn post-hoc test results with Bonferroni correction.

These objects can be used to compute effect sizes:

```r
effect_both
effect_any
```

These objects contain detailed tables of effect-size results derived from the Dunn test.

## Step 3 – Compute LOD scores for both subsets

Follow the instructions below to compute LOD scores for the two pedigree subsets:

[LOD analysis](docs/03-lod_analysis.md)

This step first computes a complete pairwise LOD score matrix and then uses it for parent retrieval at different top-*k* thresholds.

The following main objects are generated:

### 1. Pairwise LOD score matrix

```r
L
```

`L` contains the pairwise LOD scores calculated for all varieties in the dataset. Each value represents the likelihood-based score between a pair of varieties and is used as the input for the subsequent LOD parent-retrieval analysis.

The matrix is computed once and reused throughout the workflow.

### 2. LOD parent-retrieval results

```r
lod_both
lod_any
```

These objects contain the parent-retrieval results for different top-*k* values (*k* = 1, 5, 10, and 20).

To access the results for a specific top-*k* value, use, for example:

```r
lod_both$Top10
```

The returned object contains:

* detailed results;
* all_parent_scores;
* filtered_parent_scores;
* success_rate_global;
* success_rates_by_parent;
* score_stats.

## Step 4 – Perform a combined analysis on both subsets

Follow the instructions below to cperform a combined analysis for the two pedigree subsets:

[Combined analysis](docs/04-combined_analysis.md)

The following instructions generate the objects bellow:

```r
combined_both
combined_any
```

These objects contain comparison results on parent retrieval with both methods at top-10 level. In the same time: 

```r
combined_both_iqr
combined_any_iqr
```

where for Jaccard similarity metrics parent recovery was be evaluated using a similarity interval.

Each combined analysis returns a list containing

- **summary_table**

  Detailed comparison for every documented parent, including

  - Jaccard rank
  - Jaccard similarity
  - LOD rank
  - LOD score
  - recovery by each method
  - concordance status

- **global_stats**

  Summary statistics reporting

  - Both methods
  - Jaccard only
  - LOD only
  - Neither
  - Concordant
  - Discordant

as both absolute counts and percentages.

---

## Step 5 – Reproduce the manuscript figures

To generate the two main figures from the manuscript:

1. **Figure 3:** Distribution and empirical calibration of Jaccard similarity scores for documented parent–offspring pairs.
2. **Figure 4:** Distribution of Jaccard similarity scores across relationship categories in the biparental and single-parent subsets.

Follow the instructions provided in the documentation below:

[Plotting figures](docs/05-plotting.md)

Alternatively, the complete figure-generation workflow can be executed with:

```r
source("examples/reproduce_manuscript_figures.R")
```

This script runs the analyses required to generate the manuscript figures and produces the corresponding JPEG files in the repository root:

```text
fig.3-jaccard_distribution_plot.jpeg
fig.4-jaccard_structural_analysis.jpeg
```

---

## Step 6 – Reproduce the manuscript tables

To generate the four main tables reported in the manuscript:

1. **Table 1:** Retrieval performance and distribution of Jaccard similarity scores for documented parent–offspring pairs in the biparental (*n* = 266) and single-parent (*n* = 572) subsets.
2. **Table 2:** Retrieval performance and distribution of LOD scores for documented parent–offspring pairs in the biparental (*n* = 266) and single-parent (*n* = 572) subsets.
3. **Table 3:** Overlap analysis between ranking-based and structure-informed criteria for parent retrieval.
4. **Table 4:** Distribution of the number of structure-consistent candidates remaining after applying the combined LOD–Jaccard criterion to biparental offspring for which both documented parents were retained.

The documentation is provided here:

[Generating tables](docs/06-tables.md)

The complete table-generation workflow can be executed with:

```r
source("examples/reproduce_manuscript_tables.R")
```

This script runs the analyses required to generate the manuscript tables and produces the corresponding output files in the repository root: 

```text
table1_jaccard_parent_retrieval.csv
table2_lod_parent_retrieval.csv
table3_combined_lod_jaccard.csv
table4_structure_consistent_candidates.csv
```

---

# Citation

If you use this repository in your work, please cite:

Nirsha A., Spanonghe M., Rivière J., Lanterbecq D. (2026) *Improving parent identification in plant breeding using combined similarity and likelihood metrics.* (under review)

The citation will be updated once the manuscript is published.

---

# License

This repository is released under the Academic Research License.

The source code is provided to ensure scientific transparency and reproducibility of the associated publication.

Commercial use or integration into proprietary software requires prior written authorization from CARAH asbl.

For commercial licensing opportunities or other permissions, please contact: deborah.lanterbecq@condorcet.be | a.nirsha@carah.be
