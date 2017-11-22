# Create a phylogenetic tree with OTUs that have at least 2 counts throughout the dataset
## This output `newick_tree_rm2_rmN.tre` will be used for phylogenetic analysis in this manuscript

### This folder was created on May 25th, 2017 


## Prepare the fasta file from mothur and format it to run fastree 

1. `surface_PAFL_otu_pruned_RAREFIED_rm2.RData` and `OTUnames_rarefied_rm2.txt` were created with the `Create_Tree.R` script in the `code/` directory.  
2. I downloaded [bbmap](https://sourceforge.net/projects/bbmap/) and put it in my home directory.  
3. In the commandline the following was performed:  

    A. `~/bbmap/filterbyname.sh in=rm_pipe_count_16s.fa names=OTUnames_rarefied_rm2.txt out=subset_rarefied_rm2.fasta -include t`
  
    B. `sed 's/N/-/g' subset_rarefied_rm2.fasta > subset_rarefied_rm2_rmN.fasta`  

## Run Fasttree 

1. Log onto flux  
2. Make a new working directory:  For me it is `cd /scratch/lsa_fluxm/marschmi/2017_05_fasttree`  
9. Load fasttree: `module load fasttree`  
10. On flux run the `fasttree.pbs` script. It's main function is the following line of code:  


```
# Infer a tree with fasttree with the GTR+CAT 
### GTR: General time reversible model 

## Input file to fasttree = rmN_to_dash.fasta
## Output file to fasttree = newick_tree_rmN_to_dash.tre 

FastTree -gtr -nt -fastest  < subset_rarefied_rm2_rmN.fasta > newick_tree_rm2_rmN.tre
```




