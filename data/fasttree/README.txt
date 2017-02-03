# Purpose of FEb03_fasttree2 folder

This folder was created on Feb 03, 2017 by Marian Schmidt to fix fasta header names and to run fasttree for phylogenetic diversity analysis. 

# Fix header names in fasta file 
## Rename file
cp stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.rep.fasta rep_16s_seqs.fasta

## Remove the tab 

### -e in sed command allows for multiple sed commands to be run in a single sed command 

sed 's/>.*\t/>/g' rep_16s_seqs.fasta > no_tab_16s.fa 

## Remove the pipe and number at end of header name (e.g.  "|1531") 
sed 's/|.*//g' no_tab_16s.fa > rm_pipe_count_16s.fa
