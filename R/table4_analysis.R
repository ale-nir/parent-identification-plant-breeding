# 1. Extraire les deux parents documentés pour les 266 offspring
documented_parents <- varieties_both_parents %>%
  mutate(
    Parent1 = str_split_fixed(Available.Pedigree, "\\s+[Xx]\\s+", 2)[, 1],
    Parent2 = str_split_fixed(Available.Pedigree, "\\s+[Xx]\\s+", 2)[, 2]
  ) %>%
  select(Name, Parent1, Parent2)

head(documented_parents)

# 2. Pour chaque offspring, récupérer le statut des deux parents
parent_validation <- lapply(seq_len(nrow(documented_parents)), function(i) {
  
  v <- documented_parents$Name[i]
  p1 <- documented_parents$Parent1[i]
  p2 <- documented_parents$Parent2[i]
  
  lod_top10 <- lod_both$Top10$detailed_results[[v]]$top_scores
  jaccard_all <- jaccard_both$Top10$detailed_results[[v]]$all_scores
  
  candidates <- lod_top10 %>%
    left_join(jaccard_all, by = "Genotype") %>%
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
}) %>%
  bind_rows() %>%
  mutate(
    P1_Jaccard_pass = P1_Jaccard >= 0.57,
    P2_Jaccard_pass = P2_Jaccard >= 0.57,
    
    P1_Both = P1_LOD & P1_Jaccard_pass,
    P2_Both = P2_LOD & P2_Jaccard_pass,
    
    Both_parents_retained = P1_Both & P2_Both
  )
  
  table(parent_validation$Both_parents_retained)
  
table4_offspring <- table4_candidates %>%
  filter(Variety %in% parent_validation$Variety[
    parent_validation$Both_parents_retained
  ]) %>%
  group_by(Variety) %>%
  summarise(
    n_candidates = sum(Jaccard_pass, na.rm = TRUE),
    .groups = "drop"
  )

summary(table4_offspring$n_candidates)