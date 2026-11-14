
# ==========================================================
#
# CBF Gene Family Sequence Similarity Analysis
#
# This script analyses pairwise sequence distances between CBF
# genes to assess sequence identity within and between gene
# families.
#
# Breakdown of the steps:
#   - import and clean MEGA12 pairwise distance data
#   - convert evolutionary distances into percentage identity
#   - calculate descriptive statistics for each family pairing
#   - separate within-family and between-family comparisons
#
# Author: Reetta Pirttilahti
#
# ==========================================================

library(dplyr)
library(tidyr)
library(stringr)

#IMPORT DATA

distance_data <- read.csv("allfamilies_distances.csv")



#DATA FORMAT: cleaning up the data exported from MEGA12
distance_data <- distance_data %>%
  rename(seq1 = X) %>%
  mutate(
    seq1 = str_remove(seq1, "_\\(\\+1\\)$")
  )

colnames(distance_data) <- colnames(distance_data) %>%
  str_replace_all("\\.\\.[0-9]+\\.$", "") %>%
  str_replace_all("\\.\\.[0-9]+$", "") %>%
  str_replace_all("\\.", "_")

#Convert matrix to long format and filter out zeroes
dist_long <- distance_data %>%
  pivot_longer(
    cols = -seq1,
    names_to = "seq2",
    values_to = "distance"
  ) %>%
  filter(!is.na(distance))

#Extract gene family names from the gene ID in the distance matrix

family_pattern <- "IIID|IIIC|IIIB|IIIA|IAD|IVA|VB|IA|II"

dist_long <- dist_long %>%
  mutate(
    family1 = str_extract(seq1, family_pattern),
    family2 = str_extract(seq2, family_pattern)
  )

#Convert MEGA distances to identity percentage

dist_long <- dist_long %>%
  mutate(identity = (1 - distance) * 100)

#Organize family pair labels so every comparison is labelled in the same order
#for sequence similarity calculations

dist_long <- dist_long %>%
  rowwise() %>%
  mutate(
    fam_a = sort(c(family1, family2))[1],
    fam_b = sort(c(family1, family2))[2]
  ) %>%
  ungroup()



#DESCRIPTIVE STATISTICS: calculate mean distances 

family_summary <- dist_long %>%
  group_by(fam_a, fam_b) %>%
  summarize(
    mean_distance = mean(distance, na.rm = TRUE),
    mean_identity = mean(identity, na.rm = TRUE),
    min_identity = min(identity, na.rm = TRUE),
    max_identity = max(identity, na.rm = TRUE),
    n_comparisons = n(),
    .groups = "drop"
  ) %>%
  arrange(fam_a, fam_b)

#Separate within-family and between-family comparisons

within_family <- family_summary %>%
  filter(fam_a == fam_b)

between_family <- family_summary %>%
  filter(fam_a != fam_b)

