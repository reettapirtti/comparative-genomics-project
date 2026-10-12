# ==========================================================
#
# CBF Copy Number Variation (CNV) and Functional Retention
#
# This script analyses CBF gene copy number variation and
# functional retention across the A, C, and D subgenomes by:
#   - summarizing functional and non-functional gene counts
#     across species and subgenomes
#   - visualizing gene copy numbers across species and 
#     subgenomes
#   - using chi-square goodness-of-fit tests to compare gene
#     counts among subgenomes
#   - using chi-square tests of independence to assess
#     associations between functional status and subgenome
#     identity, and between functional status and gene family
#   - examining Pearson residuals to identify contributions
#     to observed differences
#
# Author: Reetta Pirttilahti
#
# ==========================================================

library(readxl)
library(tidyverse)



#   DATA IMPORT: a dataset summarizing copy numbers of functional and non-functional 
#   genes per species per sub-genome

copy_number_data <- read_excel("Copy_number_data.xlsx")



#   DATA FORMAT 


#MAIN DATASHEET

data_long <- copy_number_data %>%
  mutate(
    Subgenome = sub(".*([ACD])$", "\\1", Species_subgenome),
    Species = sub("([ACD])$", "", Species_subgenome)
  ) %>%
  select(Species, Subgenome, everything(), -Species_subgenome) %>%
  pivot_longer(
    cols = c("Functionals", "Non-functional"),
    names_to = "Status",
    values_to = "Count"
  ) %>%
  select(-`Total #`) %>%
  mutate(Count = as.numeric(Count)) %>%
  filter(!is.na(Count))


#A SUBSET FOR LATER

data_functionals <- data_long %>%
  filter(Status == "Functionals") %>%
  select(-`Gene family`) %>%
  group_by(Species, Subgenome) %>%
  mutate(`Count` = as.numeric(`Count`)) %>%
  summarise(Total = sum(`Count`, na.rm = TRUE), .groups = "drop")



#   STACKED BAR PLOT: visualizing the total number of functional genes and pseudo-
#   genes per species and per sub-genome


data_stacked <- data_long %>%
  select(Species, Subgenome, Status, Count) %>%
  group_by(Species, Subgenome, Status) %>%  
  summarise(Total = sum(Count, na.rm = TRUE), .groups = "drop")%>%
  filter(!is.na(Total) & Total > 0)

data_stacked <- data_stacked %>%
  mutate(
    Species_Subgenome = factor(
      paste(Species, Subgenome, sep = "_"),
      levels = unique(paste(Species, c("A","C","D"), sep = "_"))
    ),
    fill_color = case_when(
      Status == "Non-functional" ~ "orange",
      Subgenome == "A" ~ "darkblue",
      Subgenome == "C" ~ "darksalmon",
      Subgenome == "D" ~ "tan"
    )
  )

data_stacked$Species_Subgenome <- factor(
  data_stacked$Species_Subgenome,
  levels = c(
    "BYZ_A","BYZ_C","BYZ_D","INS_C","INS_D",
    "OW_A","OW_C","OW_D",
    "SANG_A","SANG_C","SANG_D",
    "TN1_A","TN1_C","TN1_D",
    "TN4_A","TN4_C","TN4_D"
  )
)

data_stacked$Subgenome <- factor(data_stacked$Subgenome, levels = c("A","C","D"))
species_order <- c("BYZ","INS","OW","SANG","TN1","TN4")
data_stacked$Species <- factor(data_stacked$Species, levels = species_order)


ggplot(data_stacked, aes(x = Species, y = Total, fill = fill_color)) +
  geom_bar(
    stat = "identity",
    width = 0.75,
    position = position_dodge(width = 0.8), 
    aes(group = Subgenome),
  ) +
  scale_fill_identity(
    guide = "legend",
    breaks = c("darkblue", "darksalmon", "tan"),
    labels = c("A","C","D")
  ) +
  labs(
    x = "Species",
    y = "Number of Copies",
    fill = "Subgenome"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title  = element_text(size = 11),
    panel.grid.major.x = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
  ) +
  expand_limits(y = max(data_stacked$Total) * 1.1)



#   FINDING ANSWERS IN THE DATA - 
#   STATISTICAL TESTS: 


#1) CHI-SQUARE TEST OF INDEPENDENCE


#  A) Is functional retention associated with subgenome identity (A, C, or D)?

cont_table <- data_long %>%
  group_by(Subgenome, Status) %>%
  summarise(Total = sum(Count, na.rm = TRUE)) %>%
  tidyr::pivot_wider(names_from = Status, values_from = Total, values_fill = 0)

test_functional_retention <- chisq.test(cont_table[, -1])
test_functional_retention
test_functional_retention$expected


#2) CHI SQUARE GOODNESS-OF-FIT TESTS

#  B) Are functional gene counts equal among subgenomes (A vs C vs D)?

sub_counts_functionals <- data_functionals %>% 
  group_by(Subgenome) %>%
  summarise(Total = sum(`Total`, na.rm = TRUE)) 

test_equal_counts <- chisq.test(sub_counts_functionals$Total)
test_equal_counts
test_equal_counts$residuals


#   C) Are total gene counts, including pseudogenes, distributed equally 
#   among the A, C, and D subgenomes?

sub_counts_all <- data_long %>% 
  group_by(Subgenome) %>%
  summarise(Total = sum(`Count`, na.rm = TRUE))

test_all <- chisq.test(sub_counts_all$Total)
test_all
test_all$residuals


#3) CHI SQUARE TEST OF INDEPENDENCE

#  D) Do gene families differ in the proportion of functional copies?

fam_table <- data_long %>%
  group_by(`Gene family`, Status) %>%
  summarise(Total = sum(Count, na.rm = TRUE)) %>%
  tidyr::pivot_wider(names_from = Status, values_from = Total)

test_fam <- chisq.test(fam_table[, -1])
test_fam
fam_table_results <- test_fam$residuals



#   PLOTS TO VIEW GENE FAMILY DISTRIBUTIONS


#A) All gene families to visualize the size of the families 

df_family_summary <- data_long %>%
  group_by(`Gene family`, Status) %>%
  summarise(Total = sum(Count, na.rm = TRUE), .groups = "drop") %>% 
  mutate(Status = recode(Status,
                         "Functionals" = "Functional"))

ggplot(df_family_summary,
       aes(x = `Gene family`, y = Total, fill = Status)) + 
  
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  
  geom_text(aes(label = Total),
            position = position_dodge(width = 0.75),
            vjust = -0.3,
            size = 3.5) +
  scale_fill_manual(values = c(
    "Functional" = "blue4",
    "Non-functional" = "darksalmon"
  )) +
  labs(
    x = "CBF Gene Family",
    y = "Copy Number",
    fill = "Gene Status",
    title = "Functional and non-functional copies across CBF gene families"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11,
                              margin = margin(b = 20)),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right", 
    legend.title = element_text(size= 9),
    legend.text  = element_text(size = 9),
    legend.key.size = unit(0.6, "cm")
  ) +
  expand_limits(y = max(df_family_summary$Total) * 1.1)


#B) Individual gene family plots per subgenome per species

subgenome_plots <- data_long %>% 
  filter(Status == "Functionals") %>%
  group_by(Species, Subgenome, `Gene family`) %>% 
  summarise(Total = sum(`Count`, na.rm = TRUE), .groups = "drop")

gene_families <- c("II", "IV", "IVA", "IIIC",
                   "IIIB", "IIIA", "IIID", "VB")

#The function

plot_gene_family <- function(family) {
  
  subgenome_plot <- subgenome_plots %>%
    filter(`Gene family` == family)
  
  ggplot(subgenome_plot,
         aes(x = Species, y = Total, fill = Subgenome)) +
    geom_bar(
      stat = "identity",
      width = 0.75,
      position = position_dodge(width = 0.8),
      aes(group = Subgenome)
    ) +
    scale_fill_manual(values = c(
      "A" = "#E78AC3",
      "D" = "#66C2A5",
      "C" = "#4169E1"
    )) +
    geom_text(
      aes(label = Total),
      position = position_dodge(width = 0.8),
      vjust = -0.3,
      size = 3
    ) +
    labs(
      title = paste(family, "across subgenomes"),
      x = "Species",
      y = "Total copy number"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
}

#Looping to generate each plot
subgenome_family_plots <- lapply(
  gene_families,
  plot_gene_family
)

#Naming each plot by its gene family
names(subgenome_family_plots) <- gene_families

#Displaying individual plots
subgenome_family_plots[["II"]]
subgenome_family_plots[["IV"]]
subgenome_family_plots[["IVA"]]
subgenome_family_plots[["IIIC"]]
subgenome_family_plots[["IIIB"]]
subgenome_family_plots[["IIIA"]]
subgenome_family_plots[["IIID"]]
subgenome_family_plots[["VB"]]



