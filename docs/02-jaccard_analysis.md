# Jaccard Similarity Analysis

This page describes the steps required to compute Jaccard similarity scores and perform the Jaccard-based analyses used in the manuscript.

The analysis is performed separately for the two pedigree subsets:

* `varieties_both_parents` – varieties with two documented parents;
* `varieties_any_parent` – varieties with at least one documented parent.

The workflow consists of three main analyses:

1. parent retrieval based on top-*k* Jaccard similarity;
2. comparison of Jaccard similarity among relationship groups;
3. statistical analysis of the resulting similarity distributions.

---

## 1. Compute the Jaccard similarity matrix

The pairwise Jaccard similarity matrix only needs to be computed once and is subsequently reused by all Jaccard analyses.

```r
J <- compute_jaccard_matrix_cpp(bin_matrix)
```

The resulting object `J` contains the pairwise Jaccard similarity values calculated from the binary molecular marker matrix.

---

## 2. Jaccard-based parent retrieval

The first analysis evaluates the ability of Jaccard similarity to retrieve documented parents among the most similar varieties.

Four top-*k* values are considered:

```r
top_values <- c(1, 5, 10, 20)
```

The analysis is performed independently for the two pedigree subsets:

```r
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

This generates two objects:

```r
jaccard_both
jaccard_any
```

Each object is a list containing the results for `Top1`, `Top5`, `Top10`, and `Top20`.

For example:

```r
jaccard_both$Top10
```

The results contain:

* detailed parent-retrieval results;
* parent recovery statistics;
* summary tables.

These objects are used to evaluate the performance of Jaccard similarity for parent identification at different top-*k* thresholds.

---

## 3. Structural Jaccard analysis

The second analysis evaluates whether Jaccard similarity differs according to the known biological relationship between two varieties.

Three relationship groups are considered:

* **Parent** – documented parent–offspring relationships;
* **Sibling** – varieties sharing at least one documented parent;
* **Random** – randomly selected unrelated varieties.

The relationship datasets are generated separately for each pedigree subset:

```r
jaccard_struct_both <- generate_jaccard_relationship_scores(
  data,
  J,
  varieties_both_parents
)

jaccard_struct_any <- generate_jaccard_relationship_scores(
  data,
  J,
  varieties_any_parent
)
```

This generates:

```r
jaccard_struct_both
jaccard_struct_any
```

Each object is a data frame containing one row per variety pair, together with its relationship group and Jaccard similarity score.

For example:

```r
str(jaccard_struct_both)
```

```text
'data.frame': 4854 obs. of 4 variables:
 $ Variety: chr ...
 $ Partner: chr ...
 $ Group  : chr "Parent" "Sibling" "Random" ...
 $ Jaccard: num ...
```

The resulting data frames are used as input for the statistical analysis described below.

---

## 4. Statistical analysis of Jaccard distributions

Jaccard similarity distributions are compared among the three relationship groups using a Kruskal–Wallis test followed by Dunn's post-hoc test with Bonferroni correction.

Run the analysis for each pedigree subset:

```r
results_any <- analyze_jaccard_relationships(
  jaccard_struct_any
)

results_both <- analyze_jaccard_relationships(
  jaccard_struct_both
)
```

This generates:

```r
results_both
results_any
```

Each object is a list containing:

* `summary` – descriptive statistics for each relationship group;
* `kruskal` – results of the Kruskal–Wallis test;
* `dunn` – results of Dunn's post-hoc test with Bonferroni correction.

The individual results can be accessed with:

```r
results_both$summary
results_both$kruskal
results_both$dunn
```

The same structure applies to `results_any`.

---

## 5. Calculate effect sizes

Effect sizes can be calculated from the results of Dunn's post-hoc test.

```r
effect_both <- calculate_dunn_effect_size(
  results_both$dunn,
  nrow(jaccard_struct_both)
)

effect_any <- calculate_dunn_effect_size(
  results_any$dunn,
  nrow(jaccard_struct_any)
)
```

This generates:

```r
effect_both
effect_any
```

The resulting objects contain detailed tables of effect-size estimates for the pairwise comparisons between relationship groups.

