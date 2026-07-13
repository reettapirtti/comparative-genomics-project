# ==========================================================
# 
# CBF Cluster Architecture Analysis
#
# This script analyses CBF clusters by:
#   - calculating intergenic gap sizes
#   - calculating total cluster lengths
#   - comparing sub-genomes
#   - testing statistical differences
#   - examining correlations between cluster lenght and gap size
#
# Author: Reetta Pirttilahti
# 
# ==========================================================

library(readxl)
library(dplyr)
library(tidyr)
library(grid)
library(Biostrings)
library(gridExtra)
library(grid)



#DATA IMPORT: metadata and DNA sets

cluster_data <- read_excel("CLUSTER DATA.xlsx")
clusters <- readDNAStringSet("ALL_CLUSTERS.fasta")

#DATA FORMAT
cluster_data <- cluster_data %>%
  group_by(Species, Subgenome, Gene_family) %>%
  mutate(
    Copy_no = row_number(),
    Gene_ID = paste(Species, Subgenome, Gene_family, Copy_no, sep = "_")
  ) %>%
  ungroup()

#Filter out small fragments

cluster_data <- cluster_data %>%
  select(Gene_ID, everything(), -Copy_no) %>%
  filter(Gene_family != "FRAG")



#CALCULATE GAP SIZES 

gaps <- cluster_data %>%
  filter(Gene_family != "FRAG") %>%
  arrange(Species, Subgenome, Start_bp) %>%
  group_by(Species, Subgenome) %>%
  mutate(
    Gap_start = lag(End_bp),
    Gap_end   = Start_bp,
    Gap_size  = Gap_end - Gap_start
  ) %>%
  filter(!is.na(Gap_size) & Gap_size > 0)

#Filter big gaps into separate file

big_gaps <- gaps %>%
  filter(Gap_size > 50000)

big_gaps <- big_gaps %>% 
  ungroup()

#Calculating total gap sizes in A, D, and C to compare overall sub-genome specific architecture

gaps_summary <- gaps %>%
  group_by(Subgenome) %>%
  summarise(
    total_gap_length = sum(Gap_size, na.rm = TRUE),
    n_species = n_distinct(Species),
    gap_length_per_species = total_gap_length / n_species,
    .groups = "drop"
  )

#Calculating total gaps per species per subgenome

species_gaps <- gaps %>%
  group_by(Species, Subgenome) %>%
  summarise(total_gap = sum(Gap_size, na.rm = TRUE), .groups = "drop")



#STATISTICAL TEST: do gaps differ significantly between species and between sub-genome
#Which sub-genome is driving this difference?

kruskal.test(total_gap ~ Subgenome, data = species_gaps)
pairwise.wilcox.test(species_gaps$total_gap, species_gaps$Subgenome)



#CALCULATING TOTAL CLUSTER LENGTHS 

cluster_lengths <- cluster_data %>%
  group_by(Species, Subgenome) %>%
  summarise(
    cluster_start = min(Start_bp, na.rm = TRUE),
    cluster_end = max(End_bp, na.rm = TRUE),
    cluster_length = cluster_end - cluster_start + 1,
    .groups = "drop"
  )



#CALCULATING AVG CLUSTER LENGTH PER SUB-GENOME

avg_cluster_lengths <- cluster_lengths %>%
  group_by(Subgenome) %>%
  summarise(
    avg_cluster_length = mean(cluster_length, na.rm = TRUE),
    n_species = n_distinct(Species),
    .groups = "drop"
  )



#Joining data to carry out CORRELATION TEST between cluster length and total gap size

cluster_gap_data <- cluster_lengths %>%
  left_join(species_gaps, by = c("Species", "Subgenome"))

cor.test(cluster_gap_data$cluster_length,
         cluster_gap_data$total_gap,
         method = "spearman")



#DESCRIPTIVE STATISTICS 

#Calculating variation in gaps between subgenomes 

big_gaps <- big_gaps %>% 
  group_by(Species, Subgenome) %>%
  mutate(ref_seq = paste(Species, Subgenome, sep = "_")) %>%
  ungroup()

gap_variation <- big_gaps %>%
  group_by(ref_seq) %>%
  summarise(
    mean_gap = mean(Gap_size, na.rm = TRUE),
    sd_gap   = sd(Gap_size, na.rm = TRUE),
    CV_gap   = sd_gap / mean_gap,
    n_gaps   = n(),
    .groups = "drop"
  ) %>% 
  separate(ref_seq, into = c("species", "subgenome"), sep = "_", remove = FALSE)

meta_summary <- gap_variation %>%
  group_by(subgenome) %>%   
  summarise(
    overall_mean_gap  = mean(mean_gap, na.rm = TRUE),
    overall_sd_gap    = sd(mean_gap, na.rm = TRUE),
    overall_CV_gap    = overall_sd_gap / overall_mean_gap,
    min_gap           = min(mean_gap, na.rm = TRUE),
    max_gap           = max(mean_gap, na.rm = TRUE),
    total_refseq      = n(),
    .groups = "drop"
  )

#Calculating variation in cluster length

cluster_data <- cluster_data %>%
  filter(Status == "Functional")

cluster_lengths <- cluster_data %>%
  group_by(Species, Subgenome) %>%
  summarise(
    cluster_length = max(End_bp) - min(Start_bp),
    .groups = "drop"
  )

length_variation <- cluster_lengths %>%
  group_by(Subgenome) %>%
  summarise(
    mean_length = mean(cluster_length),
    sd_length   = sd(cluster_length),
    CV_length   = sd_length / mean_length,
    .groups = "drop"
  )

#Creating a table to summarize this variation in cluster length

table_plot <- tableGrob(length_variation)

grid.newpage()
grid.draw(table_plot)
