# =============================================================================
# Reproduce manuscript tables
# =============================================================================

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

# source("setup.R")


# =============================================================================
# Table 1 – Jaccard parent-retrieval results
# =============================================================================

# Table 1 summarizes Jaccard-based parent retrieval for the biparental and
# single-parent subsets at Top-1, Top-5, Top-10 and Top-20 thresholds.
#
# Required objects:
#   jaccard_both
#   jaccard_any
#
# These objects are generated during the Jaccard analysis.
# See: docs/02-jaccard_analysis.md
# -----------------------------------------------------------------------------

top_values <- c(1, 5, 10, 20)

if (!exists("jaccard_both") || !exists("jaccard_any")) {
  stop(
    "The objects 'jaccard_both' and 'jaccard_any' are required for Table 1. ",
    "Please run the Jaccard analysis first."
  )
}


# -----------------------------------------------------------------------------
# Function to extract Table 1 statistics
# -----------------------------------------------------------------------------

extract_jaccard_table1 <- function(jaccard_results,
                                   biparental = FALSE) {

  output <- list()

  for (k in top_values) {

    result <- jaccard_results[[paste0("Top", k)]]

    # -------------------------------------------------------------------------
    # Retrieval success rates
    # -------------------------------------------------------------------------

    success_total <-
      result$success_rate_global * 100

    if (biparental) {

      success_parent1 <-
        result$success_rates_by_parent["parent1"] * 100

      success_parent2 <-
        result$success_rates_by_parent["parent2"] * 100

    }

    # -------------------------------------------------------------------------
    # Distribution of Jaccard scores
    #
    # Only documented parents successfully retained within Top-k are included.
    # -------------------------------------------------------------------------

    scores <- result$filtered_parent_scores$score

    score_stats <- c(
      Min    = min(scores, na.rm = TRUE),
      Max    = max(scores, na.rm = TRUE),
      Q1     = quantile(scores, 0.25, na.rm = TRUE),
      Q3     = quantile(scores, 0.75, na.rm = TRUE),
      Median = median(scores, na.rm = TRUE)
    )

    # -------------------------------------------------------------------------
    # Store results
    # -------------------------------------------------------------------------

    if (biparental) {

      output[[paste0("Top", k)]] <- c(
        Success_rate_total = success_total,
        Parent1_retrieval  = success_parent1,
        Parent2_retrieval  = success_parent2,
        score_stats
      )

    } else {

      output[[paste0("Top", k)]] <- c(
        Success_rate_total = success_total,
        score_stats
      )
    }
  }

  # Convert to data frame
  output <- as.data.frame(do.call(cbind, output))

  output
}


# -----------------------------------------------------------------------------
# Extract statistics for both subsets
# -----------------------------------------------------------------------------

table1_both <- extract_jaccard_table1(
  jaccard_both,
  biparental = TRUE
)

table1_any <- extract_jaccard_table1(
  jaccard_any,
  biparental = FALSE
)


# =============================================================================
# Format Table 1
# =============================================================================

# Round all numerical values to two decimal places.

table1_both <- round(table1_both, 2)
table1_any  <- round(table1_any, 2)


# Rename rows to match the manuscript

rownames(table1_both) <- c(
  "Success-rate-total (%)",
  "Parent-1-retrieval (%)",
  "Parent-2-retrieval (%)",
  "Min",
  "Max",
  "Q1",
  "Q3",
  "Median"
)

rownames(table1_any) <- c(
  "Success-rate-total (%)",
  "Min",
  "Max",
  "Q1",
  "Q3",
  "Median"
)


# =============================================================================
# Combine both subsets into one table
# =============================================================================

table1_both <- cbind(
  Metric = rownames(table1_both),
  table1_both
)

table1_any <- cbind(
  Metric = rownames(table1_any),
  table1_any
)

rownames(table1_both) <- NULL
rownames(table1_any)  <- NULL


# Add subset labels

table1_both$Subset <- "Biparental subset (n = 266)"
table1_any$Subset  <- "Single-parent subset (n = 572)"


# Reorder columns

table1_both <- table1_both[
  , c("Subset", "Metric", paste0("Top", top_values))
]

table1_any <- table1_any[
  , c("Subset", "Metric", paste0("Top", top_values))
]

# Combine the two sections

table1 <- rbind(
  table1_both,
  table1_any
)


# =============================================================================
# Export Table 1
# =============================================================================

write.csv(
  table1,
  "table1_jaccard_parent_retrieval.csv",
  row.names = FALSE
)

print(table1)

message("Table 1 was successfully generated.")

# ============================================================
# Table 2 – LOD parent retrieval
# ============================================================

extract_lod_table2 <- function(lod_results, biparental = FALSE) {

  results <- lapply(top_values, function(k) {

    result <- lod_results[[paste0("Top", k)]]

    # Overall success rate
    metrics <- tibble(
      Metric = "Success-rate-total (%)",
      Value = result$success_rate_global * 100
    )

    # Parent-specific retrieval rates for the biparental subset
    if (biparental) {
      parent_metrics <- tibble(
        Metric = c(
          "Success-rate-parent-1 (%)",
          "Success-rate-parent-2 (%)"
        ),
        Value = c(
          result$success_rates_by_parent["parent1"] * 100,
          result$success_rates_by_parent["parent2"] * 100
        )
      )

      metrics <- bind_rows(metrics, parent_metrics)
    }

    # LOD score distribution among filtered documented-parent scores
    scores <- result$filtered_parent_scores$score

    distribution <- tibble(
      Metric = c("Min", "Max", "Q1", "Q3", "Median"),
      Value = c(
        min(scores, na.rm = TRUE),
        max(scores, na.rm = TRUE),
        quantile(scores, 0.25, na.rm = TRUE),
        quantile(scores, 0.75, na.rm = TRUE),
        median(scores, na.rm = TRUE)
      )
    )

    bind_rows(metrics, distribution) %>%
      mutate(Top = paste0("Top", k))
  })

  bind_rows(results) %>%
    select(Metric, Top, Value) %>%
    tidyr::pivot_wider(
      names_from = Top,
      values_from = Value
    )
}

# Reorganize two tables

table2_both <- extract_lod_table2(
  lod_both,
  biparental = TRUE
)

table2_any <- extract_lod_table2(
  lod_any,
  biparental = FALSE
)

table2_both <- table2_both %>%
  mutate(
    Subset = "Biparental subset (n = 266)"
  ) %>%
  select(Subset, everything())

table2_any <- table2_any %>%
  mutate(
    Subset = "Single-parent subset (n = 572)"
  ) %>%
  select(Subset, everything())
  
# Combine the two sections

table2 <- bind_rows(
  table2_both,
  table2_any
) %>%
  mutate(
    across(
      starts_with("Top"),
      ~ ifelse(
        grepl("^Success-rate", Metric),
        sprintf("%.2f", .x),
        sprintf("%.0f", .x)
      )
    )
  )
  
# =============================================================================
# Export Table 2
# =============================================================================

write.csv(
  table2,
  "table2_lod_parent_retrieval.csv",
  row.names = FALSE
)

print(table2)

message("Table 2 was successfully generated.")

# ============================================================
# Table 3 – Combined LOD–Jaccard parent retrieval
# ============================================================

extract_combined_table3 <- function(combined_result, biparental = FALSE) {

  data <- combined_result$summary_table

  # In the single-parent subset, remove placeholder rows
  # corresponding to the undocumented second parent.
  if (!biparental) {
    data <- data %>%
      filter(Parent != "na")
  }

  data %>%
    mutate(
      Category = case_when(
        Selected_J & Selected_LOD ~ "Retrieved by both criteria",
        Selected_J & !Selected_LOD ~ "Retrieved by Jaccard only",
        !Selected_J & Selected_LOD ~ "Retrieved by LOD only",
        TRUE ~ "Not retrieved"
      )
    ) %>%
    count(Category, name = "Count") %>%
    mutate(
      Percentage = Count / sum(Count) * 100
    ) %>%
    select(Category, Count, Percentage)
}


# ------------------------------------------------------------
# 1. Top-10 LOD + Top-10 Jaccard
# ------------------------------------------------------------

table3_top10_both <- extract_combined_table3(
  combined_both,
  biparental = TRUE
)

table3_top10_any <- extract_combined_table3(
  combined_any,
  biparental = FALSE
)


# ------------------------------------------------------------
# 2. Top-10 LOD + Jaccard >= Q1 (0.57)
# ------------------------------------------------------------

table3_iqr_both <- extract_combined_table3(
  combined_both_iqr,
  biparental = TRUE
)

table3_iqr_any <- extract_combined_table3(
  combined_any_iqr,
  biparental = FALSE
)


# ------------------------------------------------------------
# Assemble Table 3
# ------------------------------------------------------------

category_order <- c(
  "Retrieved by both criteria",
  "Retrieved by Jaccard only",
  "Retrieved by LOD only",
  "Not retrieved"
)

table3 <- bind_rows(

  # Top-10 LOD + Top-10 Jaccard
  table3_top10_both %>%
    rename(
      Biparental_Count = Count,
      Biparental_Percentage = Percentage
    ) %>%
    left_join(
      table3_top10_any %>%
        rename(
          Single_parent_Count = Count,
          Single_parent_Percentage = Percentage
        ),
      by = "Category"
    ) %>%
    mutate(
      Criterion = "Top-10 LOD + Top-10 Jaccard"
    ),

  # Top-10 LOD + Jaccard >= Q1 (0.57)
  table3_iqr_both %>%
    rename(
      Biparental_Count = Count,
      Biparental_Percentage = Percentage
    ) %>%
    left_join(
      table3_iqr_any %>%
        rename(
          Single_parent_Count = Count,
          Single_parent_Percentage = Percentage
        ),
      by = "Category"
    ) %>%
    mutate(
      Criterion = "Top-10 LOD + Jaccard >= Q1 (0.57)"
    )
) %>%
  mutate(
    Category = factor(Category, levels = category_order)
  ) %>%
  arrange(
    factor(
      Criterion,
      levels = c(
        "Top-10 LOD + Top-10 Jaccard",
        "Top-10 LOD + Jaccard >= Q1 (0.57)"
      )
    ),
    Category
  ) %>%
  mutate(
    Biparental_Percentage = round(Biparental_Percentage, 2),
    Single_parent_Percentage = round(Single_parent_Percentage, 2)
  ) %>%
  select(
    Criterion,
    Category,
    Biparental_Count,
    Biparental_Percentage,
    Single_parent_Count,
    Single_parent_Percentage
  )


# ------------------------------------------------------------
# Export
# ------------------------------------------------------------

write.csv(
  table3,
  "table3_combined_lod_jaccard.csv",
  row.names = FALSE
)

print(table3)

message("Table 3 was successfully generated.")

# ============================================================
# Table 4 – Structure-consistent candidates after combined
# LOD–Jaccard filtering
# ============================================================

# Extract the two documented parents for the 266 biparental offspring
documented_parents <- varieties_both_parents %>%
  mutate(
    Parent1 = str_split_fixed(
      Available.Pedigree,
      "\\s+[Xx]\\s+",
      2
    )[, 1],
    Parent2 = str_split_fixed(
      Available.Pedigree,
      "\\s+[Xx]\\s+",
      2
    )[, 2]
  ) %>%
  select(Name, Parent1, Parent2)


# ------------------------------------------------------------
# Identify offspring for which both documented parents are
# retained by the combined Top-10 LOD + Jaccard >= 0.57 criterion
# ------------------------------------------------------------

parent_validation <- lapply(
  seq_len(nrow(documented_parents)),
  function(i) {

    v <- documented_parents$Name[i]
    p1 <- documented_parents$Parent1[i]
    p2 <- documented_parents$Parent2[i]

    lod_top10 <- lod_both$Top10$detailed_results[[v]]$top_scores
    jaccard_all <- jaccard_both$Top10$detailed_results[[v]]$all_scores

    candidates <- lod_top10 %>%
      left_join(
        jaccard_all,
        by = "Genotype"
      ) %>%
      mutate(
        Jaccard_pass = Score_Jaccard >= 0.57
      )

    tibble(
      Variety = v,
      Parent1 = p1,
      Parent2 = p2,

      P1_LOD = p1 %in% candidates$Genotype,
      P2_LOD = p2 %in% candidates$Genotype,

      P1_Jaccard = candidates$Score_Jaccard[
        match(p1, candidates$Genotype)
      ],
      P2_Jaccard = candidates$Score_Jaccard[
        match(p2, candidates$Genotype)
      ]
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    P1_Jaccard_pass = P1_Jaccard >= 0.57,
    P2_Jaccard_pass = P2_Jaccard >= 0.57,

    P1_Both = P1_LOD & P1_Jaccard_pass,
    P2_Both = P2_LOD & P2_Jaccard_pass,

    Both_parents_retained = P1_Both & P2_Both
  )


# ------------------------------------------------------------
# Build the candidate table used for Table 4
#
# For each biparental offspring, retain its Top-10 LOD
# candidates and determine whether each candidate also meets
# the Jaccard >= 0.57 criterion.
# ------------------------------------------------------------

table4_candidates <- lapply(
  documented_parents$Name,
  function(v) {

    lod_top10 <- lod_both$Top10$detailed_results[[v]]$top_scores
    jaccard_all <- jaccard_both$Top10$detailed_results[[v]]$all_scores

    lod_top10 %>%
      left_join(
        jaccard_all,
        by = "Genotype"
      ) %>%
      mutate(
        Variety = v,
        Jaccard_pass = Score_Jaccard >= 0.57
      )
  }
) %>%
  bind_rows()


# ------------------------------------------------------------
# Retain only offspring for which both documented parents
# were retained by the combined criterion
# ------------------------------------------------------------

table4_offspring <- table4_candidates %>%
  filter(
    Variety %in%
      parent_validation$Variety[
        parent_validation$Both_parents_retained
      ]
  ) %>%
  group_by(Variety) %>%
  summarise(
    n_candidates = sum(Jaccard_pass, na.rm = TRUE),
    .groups = "drop"
  )


# ------------------------------------------------------------
# Table 4 summary statistics
# ------------------------------------------------------------

table4 <- tibble(
  n = nrow(table4_offspring),
  Min = min(table4_offspring$n_candidates),
  Q1 = quantile(
    table4_offspring$n_candidates,
    0.25
  ),
  Median = median(table4_offspring$n_candidates),
  Mean = mean(table4_offspring$n_candidates),
  Q3 = quantile(
    table4_offspring$n_candidates,
    0.75
  ),
  Max = max(table4_offspring$n_candidates)
) %>%
  mutate(
    Q1 = as.numeric(Q1),
    Q3 = as.numeric(Q3)
  )


# Export
write.csv(
  table4,
  "table4_structure_consistent_candidates.csv",
  row.names = FALSE
)

print(table4)

message("Table 4 was successfully generated.")