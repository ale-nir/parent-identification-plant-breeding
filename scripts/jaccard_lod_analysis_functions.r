# ============================================================
# SUPPLEMENTARY CODE — Jaccard / LOD Analysis Pipeline
# ============================================================

# ------------------------------------------------------------
# 0. LIBRARIES
# ------------------------------------------------------------
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(dunn.test)
library(Rcpp)
library(readr)
library(ggsignif)

# ------------------------------------------------------------
# 1. PEDIGREE PROCESSING
# ------------------------------------------------------------

# identify varieties with both parents or at least one parent present in the dataset
analyze_pedigrees <- function(data) {

  all_names <- tolower(data$Name)

  results <- list(both_parents = character(0), any_parent = character(0))

  for (i in 1:nrow(data)) {
    current_row <- data[i, ]
    pedigree <- tolower(current_row$Available.Pedigree)

    parents <- strsplit(pedigree, " x ")[[1]]
    parents <- str_trim(parents)
    parents <- parents[parents != ""]

    p1_in_data <- !is.na(parents[1]) && parents[1] %in% all_names
    p2_in_data <- length(parents) > 1 && parents[2] %in% all_names

    if (p1_in_data && p2_in_data) {
      results$both_parents <- c(results$both_parents, current_row$Name)
    } else if (p1_in_data || p2_in_data) {
      results$any_parent <- c(results$any_parent, current_row$Name)
    }
  }

  return(list(
    both_parents = data[data$Name %in% results$both_parents, ],
    any_parent = data[data$Name %in% results$any_parent, ]
  ))
}

# ------------------------------------------------------------
# 2. JACCARD SIMILARITY
# ------------------------------------------------------------
Rcpp::sourceCpp("jaccard_cpp.cpp")
bin_matrix <- as.matrix(bin)
storage.mode(bin_matrix) <- "integer"

J <- compute_jaccard_matrix_cpp(bin_matrix)



# ------------------------------------------------------------
# 3. JACCARD PARENTAGE ANALYSIS
# ------------------------------------------------------------

analyze_parentage_jaccard <- function(data,
                                      J,
                                      varieties_with_parents,
                                      top_n) {

  ## ------------------------------------------------------------
  ## Initialisation
  ## ------------------------------------------------------------

  results <- vector("list", length(varieties_with_parents$Name))
  names(results) <- varieties_with_parents$Name

  ## noms en minuscules (une seule fois)
  data_names_lower <- tolower(data$Name)

  ## dictionnaire nom -> indice
  name_to_index <- setNames(
    seq_len(nrow(data)),
    data_names_lower
  )

  ## ------------------------------------------------------------
  ## Préparation des informations de pedigree
  ## ------------------------------------------------------------

  parent_data <- data %>%
    filter(Name %in% varieties_with_parents$Name) %>%
    mutate(
      parent1 = sapply(Available.Pedigree, function(ped) {

        parents <- strsplit(
          tolower(ped),
          " x "
        )[[1]]

        parents <- trimws(parents)

        if (length(parents) > 0)
          parents[1]
        else
          NA_character_

      }),

      parent2 = sapply(Available.Pedigree, function(ped) {

        parents <- strsplit(
          tolower(ped),
          " x "
        )[[1]]

        parents <- trimws(parents)

        if (length(parents) > 1)
          parents[2]
        else
          NA_character_

      })
    )

  ## ------------------------------------------------------------
  ## Préallocation
  ## ------------------------------------------------------------

  parent_scores_list <- vector(
    "list",
    nrow(parent_data) * 2
  )

  parent_counter <- 1

  ## ============================================================
  ## Boucle principale
  ## ============================================================

  for (i in seq_len(nrow(parent_data))) {

    var_name <- parent_data$Name[i]
    var_lower <- data_names_lower[name_to_index[[tolower(var_name)]]]

    var_id <- name_to_index[[tolower(var_name)]]

    ## ----------------------------------------------------------
    ## Jaccard déjà calculé en C++
    ## ----------------------------------------------------------

    jaccard_scores <- J[var_id, ]

    ## ----------------------------------------------------------
    ## Classement
    ## ----------------------------------------------------------

    ord <- order(
      jaccard_scores,
      decreasing = TRUE
    )

    ## retirer la variété elle-même
    ord <- ord[ord != var_id]

    ## top N
    top_idx <- ord[seq_len(top_n)]

    top_n_scores <- data.frame(
      Genotype = data$Name[top_idx],
      Score_Jaccard = jaccard_scores[top_idx],
      stringsAsFactors = FALSE
    )

    ## tous les scores
    jaccard_df <- data.frame(
      Genotype = data$Name[ord],
      Score_Jaccard = jaccard_scores[ord],
      stringsAsFactors = FALSE
    )

    current_parent_info <- parent_data[i, ]
	for (parent_num in c("parent1", "parent2")) {

      parent_name <- current_parent_info[[parent_num]]

      if (is.na(parent_name))
        next

      parent_idx <- unname(name_to_index[parent_name])

		if (is.na(parent_idx))
			next

		  if (is.null(parent_idx))
			next

      parent_score <- jaccard_scores[parent_idx]

      parent_rank <- match(parent_idx, top_idx)

      parent_in_top <- !is.na(parent_rank)

      ## stockage des scores parents
      parent_scores_list[[parent_counter]] <- data.frame(
        variety = var_name,
        parent_type = parent_num,
        parent_name = parent_name,
        in_top = parent_in_top,
        rank = parent_rank,
        score = parent_score,
        stringsAsFactors = FALSE
      )

      parent_counter <- parent_counter + 1

    }

    ## ----------------------------------------------------------
    ## Résultats détaillés
    ## ----------------------------------------------------------

    get_parent_info <- function(parent_name){

  if(is.na(parent_name)){
    return(list(
      name = NA,
      in_top = FALSE,
      rank = NA,
      score = NA
    ))
  }

  parent_idx <- unname(name_to_index[tolower(parent_name)])

  if(is.na(parent_idx)){
    return(list(
      name = parent_name,
      in_top = FALSE,
      rank = NA,
      score = NA
    ))
  }

  parent_rank <- match(parent_idx, top_idx)

  list(
    name = parent_name,
    in_top = !is.na(parent_rank),
    rank = parent_rank,
    score = jaccard_scores[parent_idx]
  )
}

    results[[var_name]] <- list(

      all_scores = jaccard_df,

      variety = var_name,

      top_scores = top_n_scores,

      parent1 = get_parent_info(current_parent_info$parent1),

      parent2 = get_parent_info(current_parent_info$parent2)

    )

  }

  ## ============================================================
  ## Fusion des parent_scores
  ## ============================================================

  parent_scores <- do.call(
    rbind,
    parent_scores_list[seq_len(parent_counter - 1)]
  )

  parent_scores_filtered <-
    parent_scores[!is.na(parent_scores$rank), ]

  ## ============================================================
  ## Statistiques
  ## ============================================================

  parent_stats <- data.frame(

    variety = names(results),

    parent1_in_top =
      sapply(results, function(x) x$parent1$in_top),

    parent2_in_top =
      sapply(results, function(x) x$parent2$in_top)

  )

  parent_stats$any_parent_in_top <-
    parent_stats$parent1_in_top |
    parent_stats$parent2_in_top

  success_rate_global <-
    mean(parent_stats$any_parent_in_top, na.rm = TRUE)

  success_rate_parent1 <-
    mean(parent_stats$parent1_in_top, na.rm = TRUE)

  success_rate_parent2 <-
    mean(parent_stats$parent2_in_top, na.rm = TRUE)

  ## ============================================================
  ## Retour
  ## ============================================================

  return(list(

    detailed_results = results,

    all_parent_scores = parent_scores,

    filtered_parent_scores = parent_scores_filtered,

    success_rate_global = success_rate_global,

    success_rates_by_parent = c(
      parent1 = success_rate_parent1,
      parent2 = success_rate_parent2
    ),

    score_stats = list(

      parent1 = summary(
        parent_scores_filtered$score[
          parent_scores_filtered$parent_type == "parent1"
        ]
      ),

      parent2 = summary(
        parent_scores_filtered$score[
          parent_scores_filtered$parent_type == "parent2"
        ]
      )

    )

  ))

}

# ------------------------------------------------------------
# 4. LOD PARENTAGE ANALYSIS
# ------------------------------------------------------------

# 4.1. LOD computation (C++)
# ------------------------------------------------------------

# LOD computation relies on a C++ implementation provided in:
# S5.2/lod_cpp.cpp
#
# To compile:
Rcpp::sourceCpp("S5.2_lod_computation/lod_cpp.cpp")

L <- compute_lod_matrix_cpp(bin_matrix)

# 4.2. Filtering of related varieties
# ------------------------------------------------------------
identify_related_varieties <- function(selected_variety, data) {

  # formatting dataset
  data <- data %>%
    mutate(
      Name = tolower(Name),
      Available.Pedigree = tolower(Available.Pedigree) %>%
        gsub("\\([^)]*\\)", "", .) %>%
        gsub("\\s+", " ", .) %>%
        str_trim()
    )

  selected_variety_clean <- tolower(selected_variety)

  # children retrieval
  children <- data %>%
    filter(str_detect(Available.Pedigree, fixed(selected_variety_clean))) %>%
    pull(Name)

  # parent retrieval
  parents_row <- data[data$Name == selected_variety_clean, ]

  parents <- if (nrow(parents_row) > 0) {
    as.character(parents_row$Available.Pedigree) %>%
      gsub("\\([^)]*\\)", "", .) %>%
      gsub("\\s+", " ", .) %>%
      str_trim()
  } else {
    ""
  }

  siblings <- character(0)

  if (length(parents) > 0 &&
      all(!is.na(parents)) &&
      all(nchar(parents) > 0) &&
      all(parents != "unknown")) {

    parent_list <- strsplit(parents, " x ")[[1]] %>%
	  str_trim()

	# Remove textual missing values
	parent_list <- parent_list[
	  !(tolower(parent_list) %in%
		  c("na", "n/a", "unknown", "?", "", "null"))
	]

    for (parent in parent_list) {

      if (nchar(parent) > 0) {

        siblings <- union(
          siblings,
          data$Name[
            str_detect(data$Available.Pedigree, fixed(parent))
          ]
        )
      }
    }

    # Remove the selected variety itself
    siblings <- setdiff(siblings, selected_variety_clean)

    # IMPORTANT:
    # A direct parent must never be classified as a sibling,
    # even if it shares one of its own parents with the offspring.
    siblings <- setdiff(siblings, parent_list)
  }

  related <- union(children, siblings)

  return(list(
    related = toupper(related),
    children = toupper(children),
    siblings = toupper(siblings)
  ))
}

# 4.3. LOD score computation
# ------------------------------------------------------------

analyze_parentage_lod <- function(data, L, varieties_list, top_k) {

  filtered_results <- list()

  parent_stats <- data.frame(
    variety = character(),
    parent1_in_top = logical(),
    parent2_in_top = logical(),
    parent1_score = numeric(),
    parent2_score = numeric(),
    stringsAsFactors = FALSE
  )

  for (var_name in varieties_list$Name) {

    var_id <- which(data$Name == var_name)
    if(length(var_id)==0) next

    var_data <- data[var_id, ]
    var_year <- var_data$Year

    pedigree <- tolower(var_data$Available.Pedigree)
    parents <- trimws(strsplit(pedigree," x ")[[1]])

    parent1 <- if(length(parents)>=1) parents[1] else NA
    parent2 <- if(length(parents)>=2) parents[2] else NA
	
	missing_values <- c("na", "n/a", "unknown", "?", "", "null")

	if (!is.na(parent1) && tolower(parent1) %in% missing_values)
	  parent1 <- NA

	if (!is.na(parent2) && tolower(parent2) %in% missing_values)
	  parent2 <- NA

    parent1_idx <- if(!is.na(parent1) && parent1 %in% tolower(data$Name))
      which(tolower(data$Name)==parent1) else NA

    parent2_idx <- if(!is.na(parent2) && parent2 %in% tolower(data$Name))
      which(tolower(data$Name)==parent2) else NA

    ## --------- only change versus original function ----------
    lod_scores <- L[var_id, ]
    ## --------------------------------------------------------

    lod_df <- data.frame(
      genotype=data$Name,
      lod_score=round(lod_scores,2),
      year=data$Year,
      stringsAsFactors=FALSE
    )

    lod_df <- lod_df[lod_df$genotype != var_name, ]
    lod_df <- lod_df[lod_df$year <= var_year, ]

    related <- identify_related_varieties(var_name,data)$related

    filtered_lod_df <- lod_df[
      !lod_df$genotype %in% related,
    ]

    filtered_lod_df <- filtered_lod_df[
      order(-filtered_lod_df$lod_score),
    ]

    filtered_topk <- head(filtered_lod_df, top_k)

    p1_in_top <- !is.na(parent1_idx) &&
      data$Name[parent1_idx] %in% filtered_topk$genotype

    p2_in_top <- !is.na(parent2_idx) &&
      data$Name[parent2_idx] %in% filtered_topk$genotype

    p1_score <- if(p1_in_top)
      filtered_topk$lod_score[
        match(data$Name[parent1_idx],filtered_topk$genotype)
      ] else NA

    p2_score <- if(p2_in_top)
      filtered_topk$lod_score[
        match(data$Name[parent2_idx],filtered_topk$genotype)
      ] else NA

    filtered_results[[var_name]] <- list(
      variety=var_name,
      parent1=if(!is.na(parent1_idx)) data$Name[parent1_idx] else NA,
      parent1_in_top=p1_in_top,
      parent1_score=p1_score,
      parent1_rank=if(p1_in_top)
        match(data$Name[parent1_idx],filtered_lod_df$genotype) else NA,
      parent2=if(!is.na(parent2_idx)) data$Name[parent2_idx] else NA,
      parent2_in_top=p2_in_top,
      parent2_score=p2_score,
      parent2_rank=if(p2_in_top)
        match(data$Name[parent2_idx],filtered_lod_df$genotype) else NA,
      topk=filtered_topk,
      excluded=c(
        related,
        data$Name[data$Year>var_year]
      )
    )

    parent_stats <- rbind(
      parent_stats,
      data.frame(
        variety=var_name,
        parent1_in_top=p1_in_top,
        parent2_in_top=p2_in_top,
        parent1_score=p1_score,
        parent2_score=p2_score,
        stringsAsFactors=FALSE
      )
    )
  }

  parent_stats$any_parent_in_top <-
    parent_stats$parent1_in_top | parent_stats$parent2_in_top

  score_stats <- list(
    parent1=list(
      n=sum(!is.na(parent_stats$parent1_score)),
      min=min(parent_stats$parent1_score,na.rm=TRUE),
      max=max(parent_stats$parent1_score,na.rm=TRUE),
      mean=mean(parent_stats$parent1_score,na.rm=TRUE),
      median=median(parent_stats$parent1_score,na.rm=TRUE),
      q1=quantile(parent_stats$parent1_score,.25,na.rm=TRUE),
      q3=quantile(parent_stats$parent1_score,.75,na.rm=TRUE),
      sd=sd(parent_stats$parent1_score,na.rm=TRUE)
    ),
    parent2=list(
      n=sum(!is.na(parent_stats$parent2_score)),
      min=min(parent_stats$parent2_score,na.rm=TRUE),
      max=max(parent_stats$parent2_score,na.rm=TRUE),
      mean=mean(parent_stats$parent2_score,na.rm=TRUE),
      median=median(parent_stats$parent2_score,na.rm=TRUE),
      q1=quantile(parent_stats$parent2_score,.25,na.rm=TRUE),
      q3=quantile(parent_stats$parent2_score,.75,na.rm=TRUE),
      sd=sd(parent_stats$parent2_score,na.rm=TRUE)
    )
  )

  list(
    detailed_results=filtered_results,
    success_rate_global=mean(parent_stats$any_parent_in_top,na.rm=TRUE),
    success_rates_by_parent=c(
      parent1=mean(parent_stats$parent1_in_top,na.rm=TRUE),
      parent2=mean(parent_stats$parent2_in_top,na.rm=TRUE)
    ),
    score_stats=score_stats
  )
}
# ------------------------------------------------------------
# 5. COMBINED ANALYSIS (Jaccard vs LOD)
# ------------------------------------------------------------

# 5.1 Jaccard top-10 with LOD top-10
# ------------------------------------------------------------
analyze_parentage_combined <- function(jaccard_results,
                                       lod_results,
                                       top_k){

  jaccard_details <- jaccard_results$detailed_results
  lod_details <- lod_results$detailed_results

  out <- vector("list", length(jaccard_details) * 2)
  k <- 1

  for(var_name in names(jaccard_details)){

    jac <- jaccard_details[[var_name]]
    lod <- lod_details[[var_name]]

    for(parent in c("parent1","parent2")){

      jp <- jac[[parent]]

      if(is.null(jp) || is.na(jp$name))
        next

      lod_rank <- lod[[paste0(parent,"_rank")]]

      out[[k]] <- data.frame(

        Variety = var_name,

        Parent = jp$name,

        Rank_J = jp$rank,

        Rank_LOD = lod_rank,

        InTop10_J =
          !is.na(jp$rank) &
          jp$rank <= top_k,

        InTop10_LOD =
          !is.na(lod_rank) &
          lod_rank <= top_k,

        stringsAsFactors = FALSE

      )

      k <- k + 1
    }
  }

  summary_df <- do.call(
    rbind,
    out[seq_len(k-1)]
  )

  summary_df$Recovered_By <- ifelse(
    summary_df$InTop10_J & summary_df$InTop10_LOD,
    "Both",
    ifelse(
      summary_df$InTop10_J,
      "Jaccard_only",
      ifelse(
        summary_df$InTop10_LOD,
        "LOD_only",
        "Neither"
      )
    )
  )

  summary_df$Discordant <-
    xor(summary_df$InTop10_J,
        summary_df$InTop10_LOD)

  stats <- list(

    Percent_J_only =
      mean(summary_df$Recovered_By=="Jaccard_only")*100,

    Percent_LOD_only =
      mean(summary_df$Recovered_By=="LOD_only")*100,

    Percent_Both =
      mean(summary_df$Recovered_By=="Both")*100,

    Percent_AtLeastOne =
      mean(summary_df$Recovered_By!="Neither")*100,

    Discordant_Cases =
      sum(summary_df$Discordant)

  )

  list(

    summary_table = summary_df,

    global_stats = stats

  )

}
# 5.2. Jaccard interquartile range with LOD top-10
# ------------------------------------------------------------

analyze_parentage_combined_modified <- function(jaccard_results,
                                                lod_results,
                                                top_k_lod,
                                                jaccard_min,
                                                jaccard_max) {

  jaccard_details <- jaccard_results$detailed_results
  lod_details <- lod_results$detailed_results

  ## ------------------------------------------------------------
  ## Summary table
  ## ------------------------------------------------------------

  out <- vector("list", length(jaccard_details) * 2)
  k <- 1

  for (var_name in names(jaccard_details)) {

    jac <- jaccard_details[[var_name]]
    lod <- lod_details[[var_name]]

    for (parent in c("parent1", "parent2")) {

      jp <- jac[[parent]]

      if (is.null(jp) || is.na(jp$name))
        next

      lod_rank <- lod[[paste0(parent, "_rank")]]

      out[[k]] <- data.frame(

        Variety = var_name,

        Parent = jp$name,

        Rank_J = jp$rank,

        Score_J = jp$score,

        Rank_LOD = lod_rank,

        InRange_J =
          !is.na(jp$score) &
          jp$score >= jaccard_min &
          jp$score <= jaccard_max,

        InTop10_LOD =
          !is.na(lod_rank) &
          lod_rank <= top_k_lod,

        stringsAsFactors = FALSE

      )

      k <- k + 1
    }
  }

  summary_df <- do.call(
    rbind,
    out[seq_len(k - 1)]
  )

  ## ------------------------------------------------------------
  ## Recovery status
  ## ------------------------------------------------------------

  summary_df$Recovered_By <- ifelse(

    summary_df$InRange_J & summary_df$InTop10_LOD,

    "Both",

    ifelse(

      summary_df$InRange_J,

      "Jaccard_only",

      ifelse(

        summary_df$InTop10_LOD,

        "LOD_only",

        "Neither"

      )
    )
  )

  summary_df$Discordant <-
    xor(summary_df$InRange_J,
        summary_df$InTop10_LOD)

  ## ------------------------------------------------------------
  ## Global statistics
  ## ------------------------------------------------------------

  stats <- list(

    Percent_J_only =
      mean(summary_df$Recovered_By == "Jaccard_only") * 100,

    Percent_LOD_only =
      mean(summary_df$Recovered_By == "LOD_only") * 100,

    Percent_Both =
      mean(summary_df$Recovered_By == "Both") * 100,

    Percent_AtLeastOne =
      mean(summary_df$Recovered_By != "Neither") * 100,

    Discordant_Cases =
      sum(summary_df$Discordant)

  )

  list(

    summary_table = summary_df,

    global_stats = stats

  )

}
# ------------------------------------------------------------
# 6. STRUCTURAL ANALYSIS
# ------------------------------------------------------------

generate_relationship_jaccard_scores <- function(data, J, varieties_list) {

  name_index <- setNames(seq_len(nrow(data)), tolower(data$Name))

  results <- vector("list", 10000)
  k <- 1

  for (var_name in varieties_list$Name) {

    cat("Processing:", var_name, "\n")

    var_id <- name_index[[tolower(var_name)]]
    if (is.null(var_id)) next

    related_info <- identify_related_varieties(var_name, data)
    children <- tolower(related_info$children)
    siblings <- tolower(related_info$siblings)

    pedigree <- tolower(data$Available.Pedigree[var_id])
    parents <- strsplit(pedigree, " x ")[[1]]
    parents <- trimws(parents)

    ## Convert textual missing values to NA
    parents[parents %in% c("", "na", "n/a", "unknown", "?", "null")] <- NA

    ## ---------------- Parents ----------------
    for (parent in na.omit(parents)) {

      parent_id <- name_index[[parent]]

      if (!is.null(parent_id)) {

        results[[k]] <- data.frame(
          Variety = var_name,
          Partner = data$Name[parent_id],
          Group = "Parent",
          Jaccard = J[var_id, parent_id],
          stringsAsFactors = FALSE
        )

        k <- k + 1
      }
    }

    ## ---------------- Siblings ----------------
    for (sibling in siblings) {

      sibling_id <- name_index[[tolower(sibling)]]

      if (!is.null(sibling_id)) {

        results[[k]] <- data.frame(
          Variety = var_name,
          Partner = data$Name[sibling_id],
          Group = "Sibling",
          Jaccard = J[var_id, sibling_id],
          stringsAsFactors = FALSE
        )

        k <- k + 1
      }
    }

    ## ---------------- Random ----------------
    excluded <- unique(tolower(c(
      var_name,
      parents,
      children,
      siblings
    )))

    available_ids <- which(!(tolower(data$Name) %in% excluded))

    if (length(available_ids) > 0) {

      random_ids <- sample(
        available_ids,
        min(5, length(available_ids))
      )

      for (rid in random_ids) {

        results[[k]] <- data.frame(
          Variety = var_name,
          Partner = data$Name[rid],
          Group = "Random",
          Jaccard = J[var_id, rid],
          stringsAsFactors = FALSE
        )

        k <- k + 1
      }
    }
  }

  do.call(rbind, results[1:(k - 1)])
}
# deep dataset analysis
analyze_jaccard_groups <- function(jaccard_df) {
  # descriptive stats
  summary_stats <- aggregate(
    Jaccard ~ Group,
    data = jaccard_df,
    FUN = function(x) {
      c(
        n = length(x),
        mean = mean(x, na.rm = TRUE),
        sd = sd(x, na.rm = TRUE),
        median = median(x, na.rm = TRUE),
        IQR = IQR(x, na.rm = TRUE)
      )
    }
  )
  summary_stats <- do.call(rbind, lapply(summary_stats, as.list))
  rownames(summary_stats) <- NULL

  # run Kruskal-Wallis computation
  kruskal_test <- kruskal.test(Jaccard ~ Group, data = jaccard_df)

  # rurn post-hoc of Dunn 
  dunn_test <- dunn.test(jaccard_df$Jaccard, jaccard_df$Group, method = "bonferroni")

  # return results
  list(
    summary = summary_stats,
    kruskal = kruskal_test,
    dunn = dunn_test
  )
}

# ------------------------------------------------------------
# 7. VISUALIZATION
# ------------------------------------------------------------

# 7.1. Plotting a dispertion plot for two datasets
# ------------------------------------------------------------

jaccard_both <- analyze_parentage_jaccard(data, J, varieties_both_parents, top_n = 30)
jaccard_any  <- analyze_parentage_jaccard(data, J, varieties_any_parent,   top_n = 30)

parent_pairs_both <- jaccard_both$all_parent_scores
parent_pairs_any <- jaccard_any$all_parent_scores


bi_pairs <- data.frame(
  pair = paste(parent_pairs_both$variety, parent_pairs_both$parent_name, sep = "--"),
  jaccard.score = parent_pairs_both$score,
  stringsAsFactors = FALSE  
)

bi_pairs <- bi_pairs[order(bi_pairs$jaccard.score), ]


uni_pairs <- data.frame(
  pair = paste(parent_pairs_any$variety, parent_pairs_any$parent_name, sep = "--"),
  jaccard.score = parent_pairs_any$score,
  stringsAsFactors = FALSE  
)

uni_pairs <- uni_pairs[order(uni_pairs$jaccard.score), ]

uni_pairs$type <- "Single-parent"
bi_pairs$type <- "Bi-parental"



df <- bind_rows(uni_pairs, bi_pairs)


# a ggplot fucntion that allows to trace a dispesion plot with an accent on "parental zone"
ggplot(df, aes(x = jaccard.score)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.02,
                 fill = "grey85",
                 color = "black") +

  geom_density(color = "black", linewidth = 0.3, adjust = 1.5) +
  annotate("rect",
           xmin = 0.61, 
           xmax = 0.69,
           ymin = 0,
           ymax = Inf,
           alpha = 0.15,
           fill = "steelblue") +
  geom_vline(xintercept = c(0.61, 0.69),
             linetype = "dashed",
             linewidth = 0.3) +
  scale_x_continuous(breaks = seq(0.35, 0.80, by = 0.05)) +

  facet_wrap(~ type, scales = "fixed", ncol = 1) +

  labs(
    x = "Jaccard similarity",
    y = "Density"
  ) +
  theme_minimal(base_size = 13)
  
# 7.2. Plotting a boxplot for structural comparison + Kluskal-Wallis integration
# ------------------------------------------------------------
# run structural analysis on both subsets
jaccard_both <- generate_relationship_jaccard_scores(data, J, varieties_both_parents)
jaccard_any  <- generate_relationship_jaccard_scores(data, J, varieties_any_parent)

#pre-process subsets so they match the needed ggplot input format
jaccard_both$Group <- as.factor(jaccard_both$Group)
jaccard_any$Group <- as.factor(jaccard_any$Group)

jaccard_both$Group <- factor(jaccard_both$Group, levels = c("Parent", "Sibling", "Random"))
jaccard_any$Group <- factor(jaccard_any$Group, levels = c("Parent", "Sibling", "Random"))

jaccard_both$Subset <- "Biparental"
jaccard_any$Subset  <- "Single-parent"

combined_data <- rbind(jaccard_both, jaccard_any)

# compute stats for both subsets
results_both <- analyze_jaccard_groups(jaccard_both)
results_any <- analyze_jaccard_groups(jaccard_any)

combined_data <- bind_rows(
  jaccard_both %>% mutate(Subset = "Biparental"),
  jaccard_any  %>% mutate(Subset = "Single-parent")
)

kw_results <- combined_data %>%
  group_by(Subset) %>%
  group_modify(~ {
    test <- kruskal.test(Jaccard ~ Group, data = .)
    tibble(
      p_value = test$p.value,
      y_pos = max(.$Jaccard) + 0.05  # facet's position
    )
  }) %>%
  ungroup() %>%
  mutate(
    label = "Kruskal–Wallis p < 2 × 10⁻¹⁶"
  )
  
  kw_results <- kw_results %>%
  mutate(x_pos = 2)

#graph dunn test + klustal-wallis

p <- ggplot(combined_data, aes(x = Group, y = Jaccard, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.15, alpha = 0.15, size = 0.6) +
  scale_fill_manual(values = c("#5DA5A4", "#C7A900", "#E67E5F")) +
  facet_wrap(~Subset, scales = "fixed") +
  labs(
    x = "Relationship category",
    y = "Jaccard similarity"
  ) +
  coord_cartesian(ylim = c(0.3, 0.88)) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(1.5, "lines"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.7)
  )

# JPEG haute résolution
ggsave(
  filename = "Figure1.jpeg",
  plot = p,
  width=29,
  height=23,
  units="cm",
  dpi = 600,
  device = "jpeg",
  bg = "white"
)

# TIFF haute résolution (recommandé pour les journaux)
ggsave(
  filename = "Figure1.tiff",
  plot = p,
   width=29,
  height=23,
  units="cm",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)
# ------------------------------------------------------------
# See S5.3 folder for usage examples
# ------------------------------------------------------------