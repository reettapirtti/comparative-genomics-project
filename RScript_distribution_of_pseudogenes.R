# ==========================================================
#
# CBF Pseudogene Distribution Analysis
#
# This script analyses the distribution of non-functional CBF
# genes within clusters by:
#   - importing a data set containing metadata on all genes 
#     in clusters, including their position and gene ID 
#   - selecting only non-functional genes
#   - calculating the position of each gene relative to its cluster,
#     allowing comparison between clusters of different lengths
#   - formatting gene metadata to prepare labels for plotting
#   - plotting pseudogene positions across subgenomes
#
# Author: Reetta Pirttilahti
#
# ==========================================================

library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)



#DATA IMPORT

cluster_data <- read_excel("cluster_data.xlsx")




#DATA FILTERING

cluster_data_NF <- cluster_data %>% 
  filter(Status == "non-functional")



#DATA FORMAT

cluster_data_NF <- cluster_data_NF %>%
  group_by(Species, Subgenome) %>%
  mutate(
    cluster_start = min(Start_bp),
    cluster_end   = max(End_bp),
    cluster_length = cluster_end - cluster_start,
    rel_position = (Start_bp - cluster_start) / cluster_length
  ) %>%
  ungroup()

cluster_data_NF <- cluster_data_NF %>%
  separate(
    Gene_ID,
    into = c("species", "subgenome", "family", "copy"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    ref_seq = paste(species, subgenome, sep = "_")) %>% 
  select(-subgenome, -species, -family, -copy)



#SET COLORS FOR PLOT

family_colors <- c("IIID"    = "blue",
                   "IIIC"    = "green",
                   "IIIA"    = "purple",
                   "IVA"     = "orange",
                   "VB"      = "magenta",
                   "IA"      = "pink",
                   "II"      = "cyan",
                   "IIIB"    = "darkgreen")



#MAKE PLOT

pseudogene_distribution <- ggplot(cluster_data_NF, 
                                  aes(x = ref_seq, y = rel_position, color = Gene_family)) +
  geom_point(size = 3, alpha = 0.8, position = position_jitter(width = 0.15)) +
  facet_wrap(~ Subgenome, scales = "free_x", nrow = 1) +  # all subgenomes in one row
  scale_color_manual(values = family_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.position = "right", 
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    title = element_text(size = 14)
  ) +
  labs(
    x = "Cluster",
    y = "Relative Position in Cluster",
    color = "Gene Family",
    title = "Pseudogene Distribution in Subgenomes"
  )

ggsave("Pseudogene_distribution.pdf", pseudogene_distribution)