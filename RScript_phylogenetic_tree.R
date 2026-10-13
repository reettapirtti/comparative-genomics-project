library(treeio)
library(ggtree)
library(readxl)
library(dplyr)



# ==========================================================
#
# Generating Phylogenetic Tree Visuals
#
# This script creates a highlighted and labelled phylogenetic tree for easy 
# visualization of previously identified and classified CBF gene clades. 
# Breakdown of the steps:
#   - import a Newick tree generated in MEGA 12 and import gene family
#     metadata 
#   - generate a fan-shaped phylogenetic tree using the R package ggtree
#   - colour terminal tips according to gene family classification 
#     in metadata data frame  
#   - highlight clades based on node numbers
#   - export the final tree as a high-resolution image
#
# Author: Reetta Pirttilahti
#
# ==========================================================



#DATA IMPORT: Newick tree export from MEGA 12 and metadata on gene families

tree <- read.tree("tree.nwk") 
meta <- read_xlsx("cluster_data - Copy.xlsx")



#DATA FORMAT

meta <- meta %>%
  mutate(
    Gene_family = as.character(ifelse(Gene_family == "IAD", "IA", Gene_family)),
    label_to_show = ifelse(Gene_family == "ref_seq", Gene_ID, NA_character_)
  )



#FINETUNING THE TREE AND HIGHLIGHTING SECTIONS BASED ON NODAL NUMBERS

p2 <- ggtree(tree, layout = "fan", open.angle = 5, branch.length = "none") %<+% meta +
  geom_tree() +
  geom_tippoint(aes(color = Gene_family), size = 2.5) +
  geom_tiplab(aes(label = label_to_show), offset = 0.9, size = 4, na.rm = TRUE) +
  scale_color_manual(
    values = c(
      "ref_seq" = "red",
      "IIID"    = "blue",
      "IIIC"    = "green",
      "IIIA"    = "purple",
      "IVA"     = "orange",
      "VB"      = "magenta",
      "IA"      = "pink",
      "II"      = "cyan",
      "IIIB"    = "darkgreen"
    ),
    name = "Gene family"  
  ) +
  theme_tree() +
  theme(
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14),
    legend.position = "bottom",
    legend.justification = c(1,0),
    legend.direction = "horizontal", 
    legend.box.margin = margin(t = 10)
  ) + 
  geom_hilight(node = 512, fill = "purple", alpha = 0.3) +
  geom_hilight(node = 300, fill = "darkgreen", alpha = 0.3) +
  geom_hilight(node = 344, fill = "pink", alpha = 0.3) +
  geom_hilight(node = 407, fill = "blue", alpha = 0.3) +
  geom_hilight(node = 396, fill = "pink", alpha = 0.3) +
  geom_hilight(node = 385, fill = "cyan", alpha = 0.3) +
  geom_hilight(node = 439, fill = "green", alpha = 0.3) +
  geom_hilight(node = 470, fill = "orange", alpha = 0.3) +
  geom_hilight(node = 487, fill = "magenta", alpha = 0.3) 

ggsave("p2_test.png", plot = p2, width = 12, height = 12, dpi = 600)


