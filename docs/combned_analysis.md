---

# Compare Jaccard and LOD results

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