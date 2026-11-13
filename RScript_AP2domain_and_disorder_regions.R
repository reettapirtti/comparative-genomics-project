
# ==========================================================
#
# AP2 Protein Domain and Disorder Analysis
#
# This script analyses AP2 protein domain architecture and
# intrinsically disordered regions across 286 CBF proteins by:
#   - using Biostrings to import, process, and export CBF protein sequences 
#   - importing protein annotation data for the CBF proteins
#     obtained using the software InterProScan
#   - identifying and removing overlapping annotations
#   - extracting AP2 domains with 15 amino acid flanking
#     regions from protein sequences
#   - using descriptive statistics to calculate mean AP2 domain
#     length, standard deviation and coefficient of variation
#   - calculating total disorder length per protein and
#     summarising disorder length across gene families
#   - using gggenes and ggplot2 to visualise protein feature
#     architecture, including AP2 domains and predicted
#     disorder regions
#   - using tidyverse packages for data manipulation and analysis

#
# Author: Reetta Pirttilahti
#
# ==========================================================


library(Biostrings)
library(tidyverse)
library(gggenes)



#DATA IMPORT 

#import protein annotation files obtained from InterProScan

first100 <- read.delim("first100.tsv", header = FALSE, sep="\t", stringsAsFactors = FALSE)
second100 <- read.delim("second100.tsv", header = FALSE, sep="\t", stringsAsFactors = FALSE)
third100 <- read.delim("third100.tsv", header = FALSE, sep="\t", stringsAsFactors = FALSE)

#import proteins

proteins <- readAAStringSet("CLUSTER_PROTEINS.fasta")



#DATA FORMAT:

#bind all the protein annotation files into one file 

protein_annotations <- bind_rows(first100, second100, third100)

#remove unnecessary columns and modify the names

protein_annotations <- protein_annotations %>%
  rename(Protein_ID = V1, Length = V3, Feature = V6, Start = V7, End = V8) %>%
  select(-V2, -V4, -V5, -V9, -V10, -V14, -V15, -V11, -V12, -V13)

protein_annotations$Protein_ID <- sub("_\\(\\+1\\)$", "", protein_annotations$Protein_ID)

#add a column for protein domain length and filter out excessively long domain annotations

protein_annotations <- protein_annotations %>%
  mutate(Length_domain = End - Start + 1)

protein_annotations <- protein_annotations %>%
  filter(Length_domain <= 200)

#create a function to remove overlapping annotations

filter_overlaps <- function(df) {
  df <- df %>% arrange(Start)
  keep <- rep(TRUE, nrow(df))
  
  for (i in seq_len(nrow(df))) {
    if (!keep[i]) next
    
    overlaps <- which(
      keep & df$Start <= df$End[i] & df$End >= df$Start[i]
    )
    
    overlaps <- setdiff(overlaps, i)
    
    if (length(overlaps) > 0) keep[overlaps] <- FALSE
  }
  
  df[keep, ]
}

#the final data frame 

clean_annotations <- protein_annotations %>%
  group_by(Protein_ID) %>%
  group_modify(~ filter_overlaps(.x)) %>%
  ungroup()



#AP2 DOMAIN EXTRACTION + FLANKING REGIONS

#add columns with metadata by extracting data from the protein_ID

clean_annotations <- clean_annotations %>%
  separate(Protein_ID, into = c("Species", "Subgenome", "Gene_family"), sep = "_", remove = FALSE)

#AP2 dataset

AP2_Dataset <- clean_annotations %>% 
  filter(Feature == "AP2 domain") 

#create a table out of CBF protein sequence data 

protein_tbl <- tibble(
  Protein_ID = sub("_\\(\\+1\\)$", "", names(proteins)),
  Sequence = as.character(proteins)
)

#add AP2 domain co-ordinates and other metadata from AP2 dataset to the CBF 
#protein sequence data so each annotation has a matching protein sequence

#include flanking regions in the extraction co-ordinates

flank_size <- 15

AP2_extract_data <- AP2_Dataset %>%
  left_join(protein_tbl, by = "Protein_ID") %>%
  mutate(
    extract_start = pmax(1, Start - flank_size),
    extract_end   = pmin(Length, End + flank_size)
  )

#extract AP2 + flanking regions 

AP2_extract_data <- AP2_extract_data %>%
  mutate(
    AP2_flank_sequence = mapply(
      function(seq, start, end) {
        substr(seq, start, end)
      },
      Sequence, extract_start, extract_end
    )
  )

#export AP2 domains and name the files according to gene family for further analysis

gene_families <- c("IIIC", "IIIB", "IA", "IIIA", "II", "IVA", "VB", "IIID")

for (fam in gene_families) {
  
  AP2_extract_family <- AP2_extract_data %>%
    filter(Gene_family == fam)
  
  AP2_flank_set <- AAStringSet(AP2_extract_family$AP2_flank_sequence)
  names(AP2_flank_set) <- AP2_extract_family$Protein_ID
  
  writeXStringSet(
    AP2_flank_set,
    paste0("AP2_flank_", fam, ".fasta")
  )
}



#DESCRIPTIVE STATISTICAL ANALYSIS 

#AVERAGE AP2 domain length 

summary_AP2 <- AP2_Dataset %>%
  summarise(
    mean_AP2_length = mean(Length_domain, na.rm = TRUE),
    sd_AP2_length   = sd(Length_domain, na.rm = TRUE),
    cv_AP2_length   = sd_AP2_length / mean_AP2_length
  )


#AVERAGE disorder length

DISORDER_dataset <- clean_annotations %>% 
  filter(Feature == "consensus disorder prediction") 

#calculate average disorder length by gene family 

summary_disorder <- DISORDER_dataset %>%
  group_by(Gene_family, Protein_ID) %>% 
  summarise(
    total_length = sum(Length_domain, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Gene_family) %>%
  summarise(mean_disorder_length = mean(total_length, na.rm = TRUE),
            sd_disorder_length   = sd(total_length, na.rm = TRUE),
            cv_disorder_length   = sd_disorder_length / mean_disorder_length
            
  )



#GRAPHS: example of visualizing AP2 domain and disorder regions in one gene family

clean_annotations <- clean_annotations %>%
  separate(
    Protein_ID,
    into = c("species", "subgenome", "family", "copy"),
    sep = "_",
    remove = FALSE
  ) 

ggenes_IIIA <- clean_annotations %>%
  filter(family == "IIIA") 

species_order <- c("INS", "OW", "BYZ", "SANG", "TN1", "TN4")

ggenes_IIIA <- ggenes_IIIA %>%
  mutate(species = factor(species, levels = species_order)) %>%
  arrange(species)

ggplot(ggenes_IIIA, aes(xmin = Start, xmax = End,
                      y = Protein_ID, fill = Feature)) +
  geom_gene_arrow(arrowhead_width = unit(2, "mm"),
                  height = 0.3, color = "black") +
  theme_minimal() +
  labs(x = "Amino acid position", y = "") +
  scale_fill_manual(values = c("AP2 domain" = "pink","AP2" = "pink", "DNA-binding domain" = "pink","consensus disorder prediction" = "gray80"))
