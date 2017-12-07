# Microhabitats shape diversity-productivity relationships in freshwater bacterial communities



###  **Authors:** Marian L. Schmidt, Bopaiah A. Biddanda, Anthony D. Weinke, Edna Chiang, Fallon Januska, Ruben Props, & Vincent J. Denef

## [Link to the Main Analysis](Final_Analysis.html)

## [Link to the Supplemental Analysis](analysis/OTU_Removal_Analysis.html)

#### Submitted to:  ADDME

  1. BioRxiv: 
  2. Journal on December XXXX, 2017  


**********

**Abstract:** Eukaryotic communities commonly display a positive relationship between biodiversity and ecosystem function (BEF). Based on current studies, it remains uncertain to what extent these findings extend to bacterial communities. An extrapolation from eukaryotic relationships would predict there to be no BEF relationships for bacterial communities because they are generally composed of an order of magnitude more taxa than the communities in most eukaryotic BEF studies. Here, we sampled surface water of a freshwater, estuarine lake to evaluate BEF relationships in bacterial communities across a natural productivity gradient. We assessed the impact of habitat heterogeneity - an important factor influencing eukaryotic BEFs - on the relationship between species richness, evenness, phylogenetic diversity, and heterotrophic productivity by sampling co-occurring free-living (more homogenous) and particle-associated (more heterogeneous) bacterial habitats. Diversity measures, and not environmental variables, were the best predictors of particles-associated heterotrophic production. There was a strong, positive, linear relationship between particle-associated bacterial richness and heterotrophic productivity that was strengthened when considering evenness. There were no observable BEF trends in free-living bacterial communities. In contrast, per-capita but not community-wide heterotrophic productivity increased across both habitats as communities were composed of taxa that were more phylogenetically clustered. This association indicates that communities with more phylogenetically related taxa have higher per-capita heterotrophic production than communities of phylogenetically distantly related taxa. Our findings show that lake heterotrophic bacterial productivity can be positively affected by evenness and richness, negatively by phylogenetic diversity, and that BEF relationships are contingent on microhabitats. These results provide a stepping stone to compare biodiversity-productivity theory developed for Eukarya to bacterial ecosystems.  

**********


## Information about this repository:  

This is the repository for the manuscript "Microhabitats shape diversity-productivity relationships in freshwater bacterial communities" written by Marian L. Schmidt, Bopaiah A. Biddanda, Anthony D. Weinke, Edna Chiang, Fallon Januska, Ruben Props, & Vincent J. Denef.  

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





This project is under the general MIT License.
