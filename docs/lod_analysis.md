# Compute the LOD matrix

Compute the LOD matrix once.

```r
L <- compute_lod_matrix_cpp(bin_matrix)
```

---

# Run the LOD analysis

For both parental subsets, compute LOD scores distribution.

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
      k
    )

  lod_any[[paste0("Top", k)]] <-
    analyze_parentage_lod(
      data,
      L,
      varieties_any_parent,
      k
    )
}
```

The returned objects contains:

- detailed parent retrieval results;
- recovery statistics;
- summary tables.

---
