# ==========================================================
#
# CBF Promoter Sequence Extraction
#
# This script extracts putative promoter sequences from CBF 
# gene clusters by:
#   - utilizing gene positional metadata together with
#     the extracted DNA sequences of CBF cluster regions
#   - extracting putative promoter sequences using the  
#     Biostrings R package
#   - dividing promoters into putative proximal and distal 
#     regions
#   - grouping promoters by gene family
#   - exporting proximal and distal sequences as FASTA files
#
# Author: Reetta Pirttilahti
#
# ==========================================================



library(Biostrings)
library(readxl)
library(dplyr)
library(gggenes)
library(IRanges)
library(writexl)



#IMPORT DATA

metadata <- read_xlsx("cluster_data.xlsx")

CLUSTERS <- readDNAStringSet("ALL_CLUSTERS.fasta")



#DATA FORMAT AND FILTERING

metadata <- metadata %>%
  mutate(
    strand = ifelse(strand == "positive", "+", "-")
  ) %>%
  mutate(ref_seq = paste(Species, Subgenome, sep = "_"))

metadata_filtered <- metadata %>%
  group_by(ref_seq) %>%
  arrange(Start_bp, .by_group = TRUE) %>%
  filter(!(row_number() == 1 & Gene_family == "IIIB")) %>%
  filter(!(row_number() == n() & Gene_family == "II")) %>%
  ungroup() %>%
  filter(Status == "Functional")



#PROMOTER SEQUENCE EXTRACTION

promoter_length <- 1500

promoters_df <- metadata_filtered %>%
  rowwise() %>%
  mutate(
    gene_seq = as.character(subseq(CLUSTERS[[ref_seq]], start = Start_bp, end = End_bp)),
    promoter_seq = if (strand == "+") {
      start_pos <- max(Start_bp - promoter_length, 1)
      as.character(subseq(CLUSTERS[[ref_seq]], start = start_pos, end = Start_bp - 1))
    } else {
      end_pos <- min(End_bp + promoter_length, length(CLUSTERS[[ref_seq]]))
      as.character(reverseComplement(subseq(CLUSTERS[[ref_seq]], start = End_bp + 1, end = end_pos)))
    }
  ) %>%
  ungroup()



#PROMOTERS ORDERED AND DIVIDED INTO DISTAL AND PROXIMAL PROMOTERS

promoters_ordered <- promoters_df %>%
  arrange(Species, Subgenome, Start_bp)

promoters_divided <- promoters_ordered %>%
  mutate(
    promoter_length = nchar(promoter_seq),
    
    proximal = substr(
      promoter_seq,
      pmax(promoter_length - 299, 1),
      promoter_length
    ),
    
    distal = substr(
      promoter_seq,
      1,
      pmax(promoter_length - 300, 1)
    )
  )



#EXPORT FAMILY SPECIFIC PROXIMAL AND DISTAL PROMOTER SEQUENCES

# Define gene families

families <- list(
  IA = c("IA", "IAD"),
  IIIA = "IIIA",
  IIIB = "IIIB",
  IIIC = "IIIC",
  IIID = "IIID",
  II = "II",
  IVA = "IVA",
  VB = "VB"
)


#A loop to extract specific sequences from the promoter files

for (family_name in names(families)) {
  
  family_data <- promoters_divided %>%
    filter(Gene_family %in% families[[family_name]])
  
  proximal_seqs <- DNAStringSet(family_data$proximal)
  names(proximal_seqs) <- family_data$Gene_ID
  
  distal_seqs <- DNAStringSet(family_data$distal)
  names(distal_seqs) <- family_data$Gene_ID
  
  writeXStringSet(
    proximal_seqs,
    paste0("proximal_promoters_", family_name, ".fasta")
  )
  
  writeXStringSet(
    distal_seqs,
    paste0("distal_promoters_", family_name, ".fasta")
  )
}

