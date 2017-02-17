In this directory, I have:

1. Subsetted the data for samples that are from the water column and the surface. (Removed MBRHP715 and MOTHJ715 because there were not enough sequences in these samples.)  
2. Rarefied the data to 2,895 sequences.  
3. Removed the sequences with 10 counts across the *entire* dataset.  
4. Wrote the data out as "nosed_merged_RAREFIED_rm10.RData".  
5. Created a file called "OTUnames_rarefied_rm10.txt", which has the name of the OTUs left in the dataset after steps 1-3.  
6. Subsetted the fasta file "rm_pipe_count_16s.fa" to only have OTUs from the subsetted data 754 OTUs by using the information in "OTUnames_rarefied_rm10.txt".  
		- Note1: "rm_pipe_count_16s.fa" was created in `data` folder (see README.md in data folder)  
		- Note2: Two commands below were run to subset the fasta file using filterbyname.sh from bbmap downloaded on February 16th, 2017  

# Run this command to create subset_rarefied_rm10.fasta
`filterbyname.sh in=rm_pipe_count_16s.fa names=OTUnames_rarefied_rm10.txt out=subset_rarefied_rm10.fasta -include t` 

# Run the following sed command to replace "N" with "-" in aligned fasta file 
sed 's/N/-/g' subset_rarefied_rm10.fasta > rmN_to_dash.fasta

7. Use the `rmN_to_dash.fasta` file to run fasttree on flux.  
8. Ran "calculate_pruned_tree.R" on flux RStudio and put output files in folder named "mpd_mntd".  
9. Read output files from "calculate_pruned_tree.R" into PrunedTree_Analysis.Rmd.  


*End*
