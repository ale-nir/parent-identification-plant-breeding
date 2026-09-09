# Parental Zone Histogram

This figure illustrates the distribution of Jaccard similarity scores for documented parent–offspring pairs in both pedigree subsets. The histogram and density curve allow visual identification of the empirical parental similarity range.

## Step 1. Run the parentage analysis

Perform the Jaccard parentage analysis for the biparental and single-parent subsets.

```r
jaccard_both <- analyze_parentage_jaccard(
    data,
    J,
    varieties_both_parents,
    top_k = 10
)

jaccard_any <- analyze_parentage_jaccard(
    data,
    J,
    varieties_any_parent,
    top_k = 10
)
```

## Step 2. Extract documented parent scores

Retrieve the Jaccard similarity scores corresponding to documented parent–offspring relationships.

```r
parent_pairs_both <- jaccard_both$all_parent_scores
parent_pairs_any  <- jaccard_any$all_parent_scores
```

## Step 3. Prepare the plotting dataset

Create a data frame containing the Jaccard scores and the pedigree subset associated with each parent–offspring pair.

```r
bi_pairs <- data.frame(
    jaccard.score = parent_pairs_both$score,
    type = "Bi-parental"
)

uni_pairs <- data.frame(
    jaccard.score = parent_pairs_any$score,
    type = "Single-parent"
)

plot_data <- dplyr::bind_rows(
    bi_pairs,
    uni_pairs
)
```

The resulting data frame contains the following columns:

| Column | Description |
|--------|-------------|
| `jaccard.score` | Jaccard similarity score of a documented parent–offspring pair. |
| `type` | Pedigree subset (`"Bi-parental"` or `"Single-parent"`). |

## Step 4. Generate the figure

Create the histogram highlighting the parental similarity zone.

```r
p <- plot_parental_zone_histogram(plot_data)
```

The function returns a **ggplot** object that can be further customized or exported using `save_publication_figure()`.

# Relationship Group Boxplot

This figure compares Jaccard similarity distributions among three biological relationship groups:

- documented parents;
- siblings;
- randomly selected unrelated varieties.

The comparison is performed separately for the biparental and single-parent pedigree subsets.

## Step 1. Generate relationship datasets

Create the relationship datasets for both pedigree subsets.

```r
jaccard_both <- generate_jaccard_relationship_scores(
    data,
    J,
    varieties_both_parents
)

jaccard_any <- generate_jaccard_relationship_scores(
    data,
    J,
    varieties_any_parent
)
```

## Step 2. Prepare the plotting dataset

Convert the relationship groups to factors and combine both pedigree subsets into a single data frame.

```r
jaccard_both$Group <- factor(
    jaccard_both$Group,
    levels = c("Parent", "Sibling", "Random")
)

jaccard_any$Group <- factor(
    jaccard_any$Group,
    levels = c("Parent", "Sibling", "Random")
)

jaccard_both$Subset <- "Biparental"
jaccard_any$Subset <- "Single-parent"

plot_data <- dplyr::bind_rows(
    jaccard_both,
    jaccard_any
)
```

The resulting data frame contains the following columns:

| Column | Description |
|--------|-------------|
| `Variety` | Focal variety. |
| `Partner` | Related or unrelated variety. |
| `Group` | Relationship category (`Parent`, `Sibling` or `Random`). |
| `Jaccard` | Jaccard similarity score. |
| `Subset` | Pedigree subset (`Biparental` or `Single-parent`). |

## Step 3. Generate the figure

```r
p <- plot_relationship_boxplot(plot_data)
```

The function returns a **ggplot** object that can be customized further or exported using `save_publication_figure()`.


# Save Publication Figures

The `save_publication_figure()` function exports a **ggplot2** object as a high-resolution image suitable for scientific publications.

Supported output formats are:

- **JPEG** (`.jpeg`)
- **TIFF** (`.tiff`, LZW compression)

By default, figures are exported using publication-quality settings:

- Resolution: **600 dpi**
- Background: **white**
- Dimensions: **29 × 23 cm**

## Usage

```r
save_publication_figure(
  plot = p,
  filename = "Figure1"
)
```

or specify the output format explicitly:

```r
save_publication_figure(
  plot = p,
  filename = "Figure1",
  format = "tiff"
)
```

## Parameters

| Parameter | Description |
|-----------|-------------|
| `plot` | A ggplot object to export. |
| `filename` | Output filename without the file extension. |
| `format` | Output format (`"jpeg"` or `"tiff"`). Default is `"jpeg"`. |
| `width` | Figure width (default: `29`). |
| `height` | Figure height (default: `23`). |
| `units` | Figure size units (default: `"cm"`). |
| `dpi` | Image resolution (default: `600`). |
