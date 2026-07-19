---

# Create pedigree subsets

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