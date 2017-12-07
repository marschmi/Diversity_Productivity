# Microhabitats shape diversity-productivity relationships in freshwater bacterial communities



###  **Authors:** Marian L. Schmidt, Bopaiah A. Biddanda, Anthony D. Weinke, Edna Chiang, Fallon Januska, Ruben Props, & Vincent J. Denef


## [Link to the Main Analysis](Final_Analysis.html)

## [Link to the Supplemental Analysis](analysis/OTU_Removal_Analysis.html)

#### Submitted to:  ADDME

  1. BioRxiv: 
  2. Journal on December XXXX, 2017  


**********

This is the repository for the manuscript "Microhabitats shape diversity-productivity relationships in freshwater bacterial communities" written by Marian L. Schmidt, Bopaiah A. Biddanda, Anthony D. Weinke, Edna Chiang, Fallon Januska, Ruben Props, & Vincent J. Denef.  The raw data can be obtained from the Sequence Read Archive at NCBI under BioProject accession number PRJNA412984.

**********


## Information about this repository:  

### **Original fastq files:**
The original Fastq files were submitted to the NCBI sequence read archive under BioProject accession number PRJNA412984.


### **File Description:**

The files in the main directory named [Final_Analysis.Rmd](https://github.com/marschmi/Diversity_Productivity/blob/master/Final_Analysis.Rmd), [Final_Analysis.html](Final_Analysis.html), and [Final_Figures/](https://github.com/marschmi/Diversity_Productivity/tree/master/Final_Figures/)` folder have all of the final analysis contents. 

Two supplementary analysis files can be found in the [analysis/](https://github.com/DenefLab/Diversity_Productivity/tree/master/analysis) folder. Specifically: Files named [OTU_Removal_Analysis.Rmd](https://github.com/marschmi/Diversity_Productivity/blob/master/analysis/OTU_Removal_Analysis.Rmd), and [OTU_Removal_Analysis.html](analysis/OTU_Removal_Analysis.html) are for the rare OTU sensitivity analysis. Files named [PrunedTree_Analysis.html](analysis/PrunedTree_Analysis.html) and [PrunedTree_Analysis.Rmd](https://github.com/DenefLab/Diversity_Productivity/blob/master/analysis/PrunedTree_Analysis.Rmd) are for different types of phylogenetic analyses. 

### **File Structure:**

```
.
├── Final_Analysis.Rmd
├── Final_Analysis.html
├── Final_Analysis.md
├── Final_Figures/
├── code/
│   ├── Create_Tree.R
│   ├── Muskegon_functions.R
│   ├── Subset_phyloseq.R
│   ├── make_otu_phyloseq.R
│   └── set_colors.R
├── data/
│   ├── PhyloTree/
│       ├── README.md
│       ├── Fasta2Phylip.pl
│       ├── OTUnames_rarefied_rm2.txt
│       ├── fasttree.pbs
│       ├── newick_tree_rm2_rmN.tre
│       ├── randomized
│       ├── rm_pipe_count_16s.fa
│       ├── subset_rarefied_rm2.fasta
│       ├── subset_rarefied_rm2_rmN.fasta
│       ├── subset_rarefied_rm2_rmN.phylip
│       └── surface_PAFL_otu_pruned_RAREFIED_rm2.RData
│   ├── fasttree/
│       ├── README.md
│       ├── fasttree.pbs
│       ├── newick_tree_16s_OTU.tre
│       ├── no_tab_16s.fa
│       ├── rep_16s_seqs.fasta
│       ├── rm_pipe_count_16s.fa
│       └── newick_tree_16s_OTU.tre
│   └── mothur/
│       ├── README.md
│       ├── mothur.1478540047.logfile
│       ├── mothur.batch.taxass
│       ├── mothur.batch.taxass.pbs
│       ├── stability.file
│       ├── stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy
│       ├── stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.rep.fasta
│       └── stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared
├── analysis/
│   ├── OTU_Removal_Analysis.Rmd
│   ├── OTU_Removal_Analysis.html
│   ├── OTU_Removal_Analysis_Figs/
│   ├── PrunedTree/
│   ├── PrunedTree_Analysis.Rmd
│   ├── PrunedTree_Analysis.html
│   └── old_analysis/
├── LICENSE
└── README.md
```

**Note:**  This project is under the general MIT License.
