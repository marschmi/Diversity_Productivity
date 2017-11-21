### This folder was created on May 25th, 2017 


1. `surface_PAFL_otu_pruned_RAREFIED_rm2.RData` and `OTUnames_rarefied_rm2.txt` were created with the `Create_Tree.R` script in the `code/` directory.  
2. I downloaded [bbmap](https://sourceforge.net/projects/bbmap/) and put it in my home directory.  
3. In the commandline the following was performed:  

    A. `~/bbmap/filterbyname.sh in=rm_pipe_count_16s.fa names=OTUnames_rarefied_rm2.txt out=subset_rarefied_rm2.fasta -include t`
  
    B. `sed 's/N/-/g' subset_rarefied_rm2.fasta > subset_rarefied_rm2_rmN.fasta`  

## Run RAxML 

4. I visited Indra Mullins Website and downloaded the `Fasta2Phylip.pl` perl script.  
5. `chmod +x Fasta2Phylip.pl`  
6. `perl Fasta2Phylip.pl subset_rarefied_rm2_rmN.fasta subset_rarefied_rm2_rmN.phylip`  
7. Log onto flux  
8. `cd /scratch/lsa_fluxm/marschmi/2017_05_RAxML`  
9. `module load raxml/8.2.8`  
10. On flux run the `raxml.pbs` script.  It's main function is the following line of code:

`raxmlHPC-MPI -s subset_rarefied_rm2_rmN.fasta -n 24samps_rm2_test -f a -m GTRGAMMA -x 777 -N 1000 -p 777 -T 20` 


Where:  

    - `raxmlHPC`: Command to run RAxML  
            - `-MPI`: Message passing interface, which allows for multithreading.  
    - `-s`: Specify the name of the alignment data file in PHYLIP or FASTA format 
    - `-n`: name of the output file  
    - `-f a`: Rapid bootstraping analysis and search for the best scoring ML tree in one program.
    - `-m`: Run the GTRGAMMA model  
    - `-x`: Specify an interger number (random seed) and turn on rapid bootstrappingan
    - `-N`: Specity the number of alternative runs on distinct starting trees. 
    - `-p`: Specify a random number seed for the parsimony inferences. this allows you to reproudce your results and will help debug.
    - `-T`: The number of cores/threads used.



to get help: `raxmlHPC -h`


### Run Fasttree 

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




