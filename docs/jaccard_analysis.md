---


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