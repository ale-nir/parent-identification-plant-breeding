# Generating manuscript tables

This page describes how to reproduce the four tables reported in the manuscript from the results generated during the Jaccard, LOD, and combined analyses.

All tables are generated automatically by the following script:

```r
source("examples/reproduce_manuscript_tables.R")
```

## Table 1 – Jaccard parent-retrieval results

Table 1 summarizes the Jaccard-based parent-retrieval results for the biparental and single-parent subsets.

Results are reported for four top-k thresholds:

- Top-1;
- Top-5;
- Top-10;
- Top-20.

For each threshold, the table reports the overall parent-retrieval success rate and the distribution of Jaccard similarity scores for the documented parent–offspring pairs.

For the biparental subset, retrieval rates are additionally reported separately for Parent 1 and Parent 2.

The following objects must have been generated during the Jaccard analysis:

```r
jaccard_both
jaccard_any
```

### Table construction

The values required for Table 1 are extracted from the summary statistics stored in these objects and assembled by top-k threshold.

The resulting table is organized into two sections:

1. Biparental subset (`n` = 266)
2. Single-parent subset (`n` = 572)

The biparental section includes overall, Parent 1, and Parent 2 retrieval rates, followed by the distribution of Jaccard similarity scores.

The single-parent section includes the overall retrieval rate and the distribution of Jaccard similarity scores.

## Table 2 – LOD parent retrieval

Table 2 summarizes the performance of LOD-based parent retrieval in the biparental and single-parent subsets.

Results are reported for four top-k thresholds:

- Top-1;
- Top-5;
- Top-10;
- Top-20.

For each threshold, the table reports the overall parent-retrieval success rate and the distribution of LOD scores for the documented parent–offspring pairs.

For the biparental subset, retrieval rates are additionally reported separately for Parent 1 and Parent 2.

The following objects must have been generated during the LOD analysis:

```r
lod_both
lod_any
```

### Table construction

The values required for Table 2 are extracted from the summary statistics stored in these objects and assembled by top-k threshold.

The resulting table is organized into two sections:

1. Biparental subset (`n` = 266)
2. Single-parent subset (`n` = 572)

The biparental section includes overall, Parent 1, and Parent 2 retrieval rates, followed by the distribution of LOD scores.

The single-parent section includes the overall retrieval rate and the distribution of LOD scores.

## Table 3 – Combined LOD–Jaccard analysis

Table 3 summarizes the overlap between LOD-based and Jaccard-based parent retrieval criteria in the biparental and single-parent subsets.

The table compares parental retrieval according to two combined criteria:

1. Top-10 LOD and Top-10 Jaccard, where a documented parent is considered retrieved when it is retained by both ranking criteria at the Top-10 threshold.
2. Top-10 LOD and Jaccard ≥ Q1 (0.57), where the LOD criterion remains restricted to the Top-10 candidates and the Jaccard criterion is replaced by a similarity threshold of 0.57.

For each criterion, parental retrieval events are classified into four mutually exclusive categories:

- Retrieved by both criteria;
- Retrieved by Jaccard only;
- Retrieved by LOD only;
- Not retrieved.

The following objects must have been generated during the combined analysis:

```r
combined_both
combined_any
combined_both_iqr
combined_any_iqr
```

### Table construction

The values required for Table 3 are extracted from the summary_table components of the combined-analysis objects.

The resulting table is organized into two sections corresponding to the two combined criteria, with biparental (`n` = 266) and single-parent (`n` = 572) subsets reported separately.

For the biparental subset, the denominator corresponds to the 532 documented parental retrieval events, representing two documented parents for each of the 266 offspring.

For the single-parent subset, the `summary_table` generated for `combined_any` and `combined_any_iqr` contains two rows per variety, including a placeholder row (`Parent = "na"`) for varieties with only one documented parent. These placeholder rows are excluded before calculating the retrieval categories, so that the analysis is based on the 572 documented parental events.

The counts and percentages are calculated from the resulting category assignments rather than directly from the `global_stats` component of the `combined_any` objects.

## Table 4 – Distribution of structure-consistent candidates

Table 4 summarizes the distribution of the number of structure-consistent candidate parents remaining after applying the combined LOD–Jaccard criterion to biparental offspring.

The criterion combines:

- the Top-10 LOD candidates; and
- a Jaccard similarity threshold of 0.57.

Only offspring for which both documented parents are retained by both criteria are included in the analysis.

The following objects must have been generated during the Jaccard and LOD analyses:

```r
lod_both
jaccard_both
varieties_both_parents
```

### Table construction

For each biparental offspring, the Top-10 LOD candidates are retrieved and matched with their corresponding Jaccard similarity scores.

Candidates with a Jaccard similarity score ≥ 0.57 are considered structure-consistent.

Only offspring for which both documented parents satisfy the combined LOD–Jaccard criterion are retained for the final analysis.

For each retained offspring, the number of candidates satisfying both criteria is counted.

The resulting distribution is summarized using:

- the number of retained offspring (`n`);
- minimum (Min);
- first quartile (Q1);
- median;
- mean;
- third quartile (Q3);
- maximum (Max).

The analysis is restricted to the biparental subset (`n` = 266).

## Output

Running:

```r
source("examples/reproduce_manuscript_tables.R")
```

generates the four manuscript tables automatically and saves the corresponding results as CSV files in the repository root:

```r
table1_jaccard_parent_retrieval.csv
table2_lod_parent_retrieval.csv
table3_combined_lod_jaccard.csv
table4_structure_consistent_candidates.csv
```