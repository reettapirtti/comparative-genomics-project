# ==========================================================
#
# Generating Consensus CBF Cluster Architectures
#
# This script generates consensus CBF cluster architectures 
# for A, C, and D subgenomes by calculating the median gene positions, 
# gene order, strand orientation, and the most 
# frequently occurring gene family at each gene position 
# across species. This is done to visualize the observed 
# structural conservation and conservation of positional organization
# within subgenomes. 
#
# Breakdown of the steps:
#   - import cluster metadata from an Excel file
#   - format gene identifiers and strand orientation
#   - filter for functional CBF genes only
#   - assign gene order based on genomic position within each cluster
#   - calculate the average number of gene positions per gene family 
#     per sub-genome across species
#   - generate consensus cluster architectures using median 
#     positions and majority gene family and strand
#   - identify large intergenic gaps between consecutive genes
#   - export the consensus data for manual plotting in BioRender
#
# Author: Reetta Pirttilahti
#
# ==========================================================


library(dplyr)
library(tidyr)
library(writexl)
library(readxl)



#DATA IMPORT (metadata file)

cluster_data <- read_excel("cluster_data.xlsx")



#DATA FORMATTING

cluster_data <- cluster_data %>%
  mutate(
    strand = ifelse(strand == "positive", "+", "-")
  )

cluster_data_NF <- cluster_data %>% 
  filter(Status == "non-functional")

cluster_data <- cluster_data %>%
  separate(
    Gene_ID,
    into = c("species", "subgenome", "family", "copy"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    ref_seq = paste(species, subgenome, sep = "_")
  ) %>%
  filter(Status == "Functional") %>%
  select(-subgenome, -species, -family, -copy)

cluster_data <- cluster_data %>%
  group_by(ref_seq) %>%
  arrange(Start_bp, .by_group = TRUE) %>%
  mutate(order_index = row_number()) %>%
  ungroup()



#CONSENSUS DATA 

avg_n_of_positions <- cluster_data %>%
  group_by(Subgenome) %>%
  summarize(avg_n_of_positions = ceiling(mean(tapply(order_index, Species, max))),
            .groups = "drop")

all_positions <- avg_n_of_positions %>%
  rowwise() %>%
  mutate(order_index = list(1:avg_n_of_positions)) %>%
  unnest(cols = c(order_index)) %>% 
  select(-avg_n_of_positions)

full_data <- all_positions %>%
  left_join(cluster_data, by = c("Subgenome", "order_index"))

consensus_cluster <- full_data %>%
  group_by(Subgenome, order_index) %>%
  summarize(
    median_start = round(median(Start_bp, na.rm = TRUE)),
    median_end   = round(median(End_bp, na.rm = TRUE)),
    majority_strand = names(sort(table(strand), decreasing = TRUE))[1],
    majority_family = names(sort(table(Gene_family), decreasing = TRUE))[1],
    majority_count = max(table(Gene_family)),
    .groups = "drop"
  )


#Add a column which recognizes "big gaps" for ease of plotting

consensus_cluster <- consensus_cluster %>%
  group_by(Subgenome) %>%
  arrange(order_index, .by_group = TRUE) %>%
  mutate(
    gap_size = median_start - lag(median_end),
    Gaps = ifelse(gap_size > 50000, gap_size, NA)
  ) %>%
  ungroup()


#Save file for manual plotting in BioRender 

write_xlsx(consensus_cluster, "Consensus cluster.xlsx")
