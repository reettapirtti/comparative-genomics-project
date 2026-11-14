# ==========================================================
#
# CBF Promoter Motif Visualisation
#
# This script generates visuals for predicted transcription  
# factor binding sites in selected CBF gene promoter sequences by:
#   - importing FIMO-predicted binding site data
#   - combining motif data across gene families
#   - selecting promoter sequences of interest
#   - preparing motif data for visualisation
#   - plotting motif positions along promoter sequences using ggplot
#
# Author: Reetta Pirttilahti
#
# ==========================================================



library(tidyverse)
library(readxl)



#DATA IMPORT: import FIMO-predicted transcription factor binding sites
#for each gene family, identified using a database from Plant TFDB

motifs_IIIA <- read_xlsx("IIIA_FIMO_clean_filtered5.xlsx")
motifs_IA <- read_xlsx("IA_FIMO_clean_filtered5.xlsx")
motifs_IIIC <- read_xlsx("IIIC_FIMO_clean_filtered5.xlsx")
motifs_IVA <- read_xlsx("IVA_FIMO_clean_filtered5.xlsx")
motifs_VB <- read_xlsx("VB_FIMO_clean_filtered5.xlsx")



#DATA FORMAT: combine motif data and extract motif family names

combined_motifs <- bind_rows(
  motifs_IIIA,
  motifs_IA,
  motifs_IIIC,
  motifs_IVA,
  motifs_VB
)

combined_motifs <- combined_motifs %>%
  mutate(motif = str_extract(motif_alt_id, "[^.]+$"))



#SELECT GENES WITH INTERESTING MOTIF PATTERNS FOR VISUALIZATION

seqs_IA <- c("TN4_A_IA_1",
             "TN4_A_IA_2",
             "TN4_A_IA_3",
             "TN4_A_IA_4",
             "TN4_A_IA_5"
)

seqs_IVA <- c("TN4_C_IVA_3",
              "TN4_C_IIIC_1",
              "TN4_C_VB_1"
) 



#PREPARE MOTIF DATA

plot_IA <- combined_motifs %>%
  filter(sequence_name %in% seqs_IA) %>%
  mutate(
    start = as.numeric(start),
    stop = as.numeric(stop),
    sequence_name = factor(sequence_name, levels = rev(seqs_IA))
  )

plot_IVA <- combined_motifs %>%
  filter(sequence_name %in% seqs_IVA) %>%
  mutate(
    start = as.numeric(start),
    stop = as.numeric(stop),
    sequence_name = factor(sequence_name, levels = rev(seqs_IVA))
  )

#Rename motifs, e.g. ERF 1-9 are all called ERF for simplicity

plot_IVA$motif <- gsub("^ERF[0-9]+", "ERF", plot_IVA$motif)
plot_IVA$motif[plot_IVA$motif %in% c("MYB_RELATED", "MYB_related")] <- "MYB-related"



#GENERATE DATA FRAMES FOR PLOTTING 

df_IA <- tibble(
  sequence_name = factor(rev(seqs_IA), levels = rev(seqs_IA)),
  xmin = 1,
  xmax = 1300,
  y = factor(rev(seqs_IA), levels = rev(seqs_IA))
)

df_IVA <- tibble(
  sequence_name = factor(rev(seqs_IVA), levels = rev(seqs_IVA)),
  xmin = 1,
  xmax = 1300,
  y = factor(rev(seqs_IVA), levels = rev(seqs_IVA))
)



#GENERATE VISUALS USING GGPLOT 

ggplot() +
  geom_segment(
    data = df_IA,
    aes(x = xmin, xend = xmax, y = y, yend = y),
    linewidth = 1.2,
    colour = "grey40"
  ) +
  geom_rect(
    data = plot_IA,
    aes(
      xmin = start,
      xmax = stop,
      ymin = as.numeric(sequence_name) - 0.2,
      ymax = as.numeric(sequence_name) + 0.2,
      fill = motif
    )
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = "Promoter position (bp)",
    y = NULL,
    fill = "Motif family"
  ) +
  scale_x_continuous(
    limits = c(1, 1300),
    breaks = seq(0, 1300, by = 200),
    expand = c(0.01, 0.01)
  ) +
  scale_y_discrete() +
  scale_fill_manual(values = c(
    "chartreuse","blue","lightsalmon","orange","magenta","forestgreen",
    "purple","firebrick1","honeydew4","red4","cadetblue1","yellow2",
    "black","cyan4"
  )) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right"
  )

ggplot() +
  geom_segment(
    data = df_IVA,
    aes(x = xmin, xend = xmax, y = y, yend = y),
    linewidth = 1.2,
    colour = "grey40"
  ) +
  geom_rect(
    data = plot_IVA,
    aes(
      xmin = start,
      xmax = stop,
      ymin = as.numeric(sequence_name) - 0.2,
      ymax = as.numeric(sequence_name) + 0.2,
      fill = motif
    )
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = "Promoter position (bp)",
    y = NULL,
    fill = "Motif family"
  ) +
  scale_x_continuous(
    limits = c(1, 1300),
    breaks = seq(0, 1300, by = 200),
    expand = c(0.01, 0.01)
  ) +
  scale_y_discrete() +
  scale_fill_manual(values = c(
    "chartreuse","blue","magenta","cadetblue1","black","yellow",
    "forestgreen"
  )) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right"
  )

