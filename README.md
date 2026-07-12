# comparative-genomics-project
# general overview
Independent comparative genomics analysis of CBF gene clusters in oat, developed as part of my Master's thesis: "Structural Evolution and Subgenome-Specific Divergence of CBF Gene Clusters in Wild and Cultivated Oat". I developed bioinformatics pipelines and performed statistical analyses of genomic datasets using RStudio.

---

# overview of the data I used 
Around 1Mb size genomic regions containing CBF clusters from chromosome 5 of the A, C, and D subgenomes from 6 different species of oat: the tetraploid _Avena insularis_ (only containing the C and the D subgenomes), an old winter _Avena sativa_ cultivar, the Mediterranean _Avena byzantina_, a spring _Avena sativa_ cultivar cv. Sang, and two wild Turkish _Avena sterilis_: _Avena sterilis_ sp. TN4, and _A. sterilis_ sp. TN1. 

The genomic regions were annotated for CBF genes in CLC Main Workbench for a total of 298 CBF genes with intact reading frames, and data was collected on gene orientation, gene length, and gene location in the clusters. Genes were grouped together into clades through iterative phylogenetic reconstruction and the gene families were named in accordance with reference sequences in the literature.  

---

## repository-contents

*consensus gene clusters: 

The positional organization of gene families was found to be highly conserved within subgenomes across species. In order to visualize the structural conservation and variation between A, C, and D subgenome clusters, consensus values summarizing gene order, positions, strand orientation, and gene family were generated. The consensus data was used to manually plot genes on a vertical, scaled axis in BioRender to generate nice figures. 

*gaps and intergenic variation, statistical analysis and data handling: 

While positional organization of genes was highly conserved, cluster length varied significantly both between species and between sub-genomes due to intergenic gap variation. To quantify and interpret these differences, data on intergenic gap length was collected and handled in RStudio. 

*copy number variation and functional retention: 

*distribution of pseudogenes: 

*promoter region extraction: 

##packages
*dplyr
*tidyr
*writexl
*readxl
*tidyverse 
*BiocManager
*GenomicRanges 
*Biostrings 
