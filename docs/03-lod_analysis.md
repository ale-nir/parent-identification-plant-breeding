# Compute the LOD matrix

Compute the LOD matrix once.

```r
L <- compute_lod_matrix_cpp(bin_matrix)
```

---

# Run the LOD analysis

For both parental subsets, compute LOD scores distribution.

If `top_values` has not already been created during the Jaccard analysis, define it as follows:

```r
top_values <- c(1, 5, 10, 20)
```

Then: 

```r
lod_both <- list()
lod_any  <- list()

for(k in top_values){

  lod_both[[paste0("Top", k)]] <-
    analyze_parentage_lod(
      data,
      L,
      varieties_both_parents,
      top_n = k
    )

  lod_any[[paste0("Top", k)]] <-
    analyze_parentage_lod(
      data,
      L,
      varieties_any_parent,
      top_n = k
    )
}
```

This generates two objects:

```r
lod_both
lod_any
```

Each object is a list containing the results for `Top1`, `Top5`, `Top10`, and `Top20`.

For example:

```r
lod_both$Top10
```

The results contain:

* detailed results;
* all_parent_scores;
* filtered_parent_scores;
* success_rate_global;
* success_rates_by_parent;
* score_stats.

---
