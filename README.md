# comparative-genomics-project
## general-overview
Independent comparative genomics analysis of CBF gene clusters in oat, developed as part of my Master's thesis: "Structural Evolution and Subgenome-Specific Divergence of CBF Gene Clusters in Wild and Cultivated Oat". I developed bioinformatics pipelines and performed statistical analyses of genomic datasets using RStudio.

This repository is a portfolio of my skills in comparative genomic analysis.

---

## overview-of-the-data 
Around 1Mb size genomic regions containing CBF clusters from chromosome 5 of the A, C, and D subgenomes from 6 different species of oat: the tetraploid _Avena insularis_ (only containing the C and the D subgenomes), an old winter _Avena sativa_ cultivar, the Mediterranean _Avena byzantina_, a spring _Avena sativa_ cultivar cv. Sang, and two wild Turkish _Avena sterilis_: _Avena sterilis_ sp. TN4, and _A. sterilis_ sp. TN1. 

The genomic regions were annotated for CBF genes in CLC Main Workbench for a total of 298 CBF genes with intact reading frames, and metadata on gene orientation, gene length, and gene location in the clusters was collected. Genes were grouped together into clades through iterative phylogenetic reconstruction and the gene families were named in accordance with reference sequences in the literature.  

The data has not been included in this repository for confidentiality reasons.

---

## repository-contents

### SEQUENCE ANALYSIS AND PHYLOGENETICS

#### ~phylogenetic-tree

To exhibit the size of gene families and their relationships to known CBF genes, a neighbor-joining tree using translated full-length CBF gene sequences was created in MEGA12 and imported into R as a Newick file. The visuals for the topology-only phylogenetic tree were created using the package ggtree. 

#### ~distances

Translated CBF genes were aligned using ClustalW and pairwise distances were calculated in MEGA12. The results were processed in R to to assess sequence identity within and between gene families. CBF protein sequences within the same clade showed high sequence similarity (mean ≤82%), while genes in different clades exhibited lower identities ranging from 41-68%.

#### ~AP2-domain-and-disorder-regions

The predicted functionality of identified full-length CBF genes was assessed by the presence of a full-length AP2 domain. AP2 domain regions were extracted from 298 translated CBF genes by utilizing protein annotation data obtained using the software InterProScan and by utilizing the Rpackage Biostrings. Descriptive statistics were used to summarise AP2 domain lengths across the identified CBF proteins, including the mean, standard deviation, and coefficient of variation. 


### STRUCTURAL ANALYSIS

#### ~consensus-gene-clusters: 

The positional organization of gene families was found to be highly conserved within subgenomes across species. In order to visualize the observed structural conservation and variation between A, C, and D subgenome clusters, consensus values summarizing gene order, positions, strand orientation, and gene family were generated. The consensus data was used to manually plot genes on a vertical, scaled axis in BioRender to generate nice figures. 

#### ~copy-number-variation-(CNV)-and-functional-retention: 

To determine whether gene family expansions or losses were biased among subgenomes or species, copy number variation and functional retention were quantified. The dataset for this summarized the number of functional and non-functional genes from each gene family per species per subgenome. Non-functional genes, or pseudogenes, were defined as sequences that were highly similar to their functional counterparts (<80% identity, and over <550bp in length), but their CDS was not intact. 

#### ~gaps-and-intergenic-variation, statistical-analysis-and-data-handling: 

While positional organization of genes was highly conserved, cluster length varied significantly both between species and between sub-genomes due to intergenic gap variation. To quantify and interpret these differences, metadata collected on gene position, orientation, and length, was used in RStudio. Statistical analyses: Kruskal–Wallis test, pairwise Wilcoxon rank-sum tests (post hoc), Spearman's rank correlation, and descriptive statistics (mean, standard deviation, coefficient of variation, minima/maxima, totals).

#### ~distribution-of-pseudogenes

The distribution of CBF pseudogenes across the CBF clusters was observed to be non-random. To visualise this pattern and compare pseudogene locations across clusters of varying lengths, gene positions were normalised relative to the length of their respective clusters. The relative positions of non-functional CBF genes were then plotted across the subgenomes, with gene families distinguished by colour.  

### PROMOTER ANALYSIS 

#### ~promoter-region-extraction 

In RStudio, 1500bp of promoter sequence upstream of the start codon was extracted for each gene. The R package Biostrings was used for batch sequence manipulation and extraction. The 1500bp extracted region was separated into a putative core region for core promoter analysis (the last 300bp of the extraction) and a proximal/distal promoter
region to search for binding sites for regulatory TFs (the first 1200bp). 

#### ~promoter-motif-visualization

Putative transcription factor binding site motifs within the extracted distal promoter regions were identified using FIMO, with motif databases obtained from JASPAR and PlantTFDB. The results were downloaded as TSV files and batch processed in Rstudio using the Rpackages Biostrings and IRanges. FIMO results were filtered by removing overlapping motif matches, retaining the most statistically significant hit (lowest p-value) and the longest motif length. Promoter motif schematics were generated in RStudio using the Rpackage gggenes.
