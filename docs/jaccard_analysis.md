# Compute the Jaccard similarity matrix

The similarity matrix only needs to be computed once.

```r
J <- compute_jaccard_matrix_cpp(bin_matrix)
```

---

# Run the Jaccard analysis

For both parental subsets, compute Jaccard scores distribution.

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

# Structural Jaccard Analysis

This analysis compares Jaccard similarity scores among biologically related and unrelated varieties. Three relationship groups are considered:

- documented parents;
- siblings;
- randomly selected unrelated varieties.

Generate the relationship dataset for each pedigree subset.

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

Each function returns a data frame containing one row per variety pair.

Example:

```r
str(jaccard_struct_both)
```

```
'data.frame': 4854 obs. of 4 variables:
 $ Variety: chr ...
 $ Partner: chr ...
 $ Group  : chr "Parent" "Sibling" "Random" ...
 $ Jaccard: num ...
```

The relationship groups are:

- **Parent** – documented parent–offspring relationships;
- **Sibling** – varieties sharing at least one documented parent;
- **Random** – randomly selected unrelated varieties.

---

# Statistical Analysis

Compare Jaccard similarity distributions among relationship groups.

```r
results_any <- analyze_jaccard_relationships(
  jaccard_struct_any
)

results_both <- analyze_jaccard_relationships(
  jaccard_struct_both
)
```

Each function returns a list containing:

- `summary` – descriptive statistics for each relationship group;
- `kruskal` – Kruskal–Wallis test;
- `dunn` – Dunn post-hoc test with Bonferroni correction.

Example:

```r
results_both$summary
results_both$kruskal
results_both$dunn
```