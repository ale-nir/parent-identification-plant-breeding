# parent-identification-plant-breeding
<i>This repository accompanies the manuscript which is actually under review</i>

"Improving parent identification in plant breeding using combined similarity and likelihood metrics"
Alexandra Nirsha¹, John Riviere², Martin Spanoghe², Deborah Lanterbecq¹,²,³

¹ Centre pour l’Agronomie et l’Agro-Industrie de la Province de Hainaut (asbl CARAH), 7800 Ath, Belgium

² Laboratoire de Biotechnologie et Biologie Appliquée, Haute Ecole Provinciale de Hainaut-CONDORCET, 7000 Mons, Belgium

³ Hainaut Analyses (HA), 7000 Mons, Belgium

# The repo content:

- implementation of the Jaccard similarity metric
- implementation of the binary LOD score
- scripts reproducing the figures of the manuscript
- anonymized example dataset

---

# Parentage Analysis Pipeline

This package provides a complete workflow for evaluating parent–offspring relationships from binary molecular marker data using Jaccard similarity, LOD scores, and their combined analysis.

The typical workflow consists of:

1. Loading the project and required functions.
2. Creating pedigree subsets.
3. Computing similarity matrices.
4. Running Jaccard analyses.
5. Running LOD analyses.
6. Comparing both methods using the combined analysis.

---

# Installation

## Clone the repository

```bash
git clone https://github.com/ale-nir/parent-identification-plant-breeding.git
```
or download the repository as a ZIP archive from GitHub and extract it.

---

## Open the project

Open the project in RStudio (recommended) or set the working directory to the repository root.

```r
setwd("path/to/repository")
```

---

## Load all package functions

Source the setup script.

```r
source("R/setup.R")
```

This script automatically loads every function required by the analysis pipeline as well as datasets.

---

# Step 1 – Create pedigree subsets

Create the pedigree subsets used during validation.

```r
ped <- create_pedigree_subsets(data)

varieties_both_parents <- ped$biparental_subset
varieties_any_parent  <- ped$uniparental_subset
```

Two validation datasets are produced:

- **Biparental subset**: varieties with two documented parents.
- **Uniparental subset**: varieties with at least one documented parent.

---

# Step 2 – Compute the Jaccard similarity matrix

The similarity matrix only needs to be computed once.

```r
J <- compute_jaccard_matrix_cpp(bin_matrix)
```

---

# Step 3 – Run the Jaccard analysis

Run the analysis for each pedigree subset.

```r
top_values <- c(1, 5, 10, 20)

jaccard_both <- list()
jaccard_any  <- list()

for(k in top_values){

  jaccard_both[[paste0("Top", k)]] <-
    analyze_parentage_jaccard(
      data,
      J,
      varieties_both_parents,
      top_k = k
    )

  jaccard_any[[paste0("Top", k)]] <-
    analyze_parentage_jaccard(
      data,
      J,
      varieties_any_parent,
      top_k = k
    )

}
```

The returned objects contains:

- detailed parent retrieval results;
- recovery statistics;
- summary tables.

---

# Step 4 – Compute the LOD matrix

Compute the LOD matrix once.

```r
L <- compute_lod_matrix_cpp(bin_matrix)
```

---

# Step 5 – Run the LOD analysis

```r
top_values <- c(1, 5, 10, 20)

lod_both <- list()
lod_any  <- list()

for(k in top_values){

  lod_both[[paste0("Top", k)]] <-
    analyze_parentage_lod(
      data,
      L,
      varieties_both_parents,
      top_k = k
    )

  lod_any[[paste0("Top", k)]] <-
    analyze_parentage_lod(
      data,
      L,
      varieties_any_parent,
      top_k = k
    )

}
```

The returned objects contains:

- detailed parent retrieval results;
- recovery statistics;
- summary tables.

---

# Step 6 – Compare Jaccard and LOD results

## Top-k comparison

Compare parent recoveries using the same Top-k threshold.

```r
combined_both <- analyze_parentage_combined(
    jaccard_both$Top10,
    lod_both$Top10,
    jaccard_selection = "top",
    top_k = 10
)

combined_any <- analyze_parentage_combined(
    jaccard_any$Top10,
    lod_any$Top10,
    jaccard_selection = "top",
    top_k = 10
)
```

---

## Jaccard score interval comparison

Instead of selecting the Top-k Jaccard candidates, parent recovery can be evaluated using a similarity interval.

For example, using the empirical first quartile (0.57):

```r
combined_both_iqr <- analyze_parentage_combined(
    jaccard_both$Top10,
    lod_both$Top10,
    jaccard_selection = "range",
    top_k = 10,
    jaccard_min = 0.57,
    jaccard_max = 1
)

combined_any_iqr <- analyze_parentage_combined(
    jaccard_any$Top10,
    lod_any$Top10,
    jaccard_selection = "range",
    top_k = 10,
    jaccard_min = 0.57,
    jaccard_max = 1
)
```

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
