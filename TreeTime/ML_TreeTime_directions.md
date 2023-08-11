Building a maximum-likelihood TimeTree using RAxML output
---

This tutorial outlines methods for building a maximum-likelihood time-scaled phylogenies described [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5758920/) with tutorials listed [here](https://github.com/neherlab/treetime). I will use the same tree files as done in building the ML trees: picornavirales full genomes, picornavirales p1 region, picornavirales rdrp region, picornaviridae full genomes, picornaviridae p1 region, picornaviridae rdrp region. Each is in their own folder in the master folder [here](). I chose this method over BEAST as my datasets are huge. 

---
General instructions: 
___
1. I took the FBP support files from the RAxML ML-tree outputs and opened them in FigTree in order to get the newick files. Then, using the metedata csv files used in the RAxML ML-tree folders, I created separate ones that only had the accession number and the full date of collection in yyyy-mm-dd format. 
---

2. Next, I added the cleaned and aligned fasta files that were used in the original ML-tree analyses

---
3. Finally, after installing following the instructions in the GitHub linked above, I used the command line to cd into each folder with the newick file, the aligned fasta file, and the metadata file with the accessions and full dates, I ran the following code for each tree type, with a general TimeTree,and one with a coalescent model (which I would have done with a Bayesian tree based on literature)

picornavirales all
```
cd Desktop/mada-bat-picornavirus/TreeTime/picornavirales_all

treetime --aln picornavirales_all_align.fasta --tree picornavirales_all_rooted --dates picornavirales_all_treetime_date.csv --outdir picornavirales_all_rooted_out

treetime --aln picornavirales_all_align.fasta --tree picornavirales_all_rooted --dates picornavirales_all_treetime_date.csv --outdir picornavirales_all_rooted_skyline_out --coalescent skyline

treetime --aln picornavirales_all_align.fasta --tree picornavirales_all_unrooted --dates picornavirales_all_treetime_date.csv --outdir picornavirales_all_unrooted_out

treetime --aln picornavirales_all_align.fasta --tree picornavirales_all_unrooted --dates picornavirales_all_treetime_date.csv --outdir picornavirales_all_unrooted_skyline_out --coalescent skyline

```


picornavirales full
```
cd Desktop/mada-bat-picornavirus/TreeTime/picornavirales_full

treetime --aln picornavirales_full_align.fasta --tree picornavirales_full_rooted --dates picornavirales_full_treetime_date.csv --outdir picornavirales_full_rooted_out

treetime --aln picornavirales_full_align.fasta --tree picornavirales_full_rooted --dates picornavirales_full_treetime_date.csv --outdir picornavirales_full_rooted_skyline_out --coalescent skyline

treetime --aln picornavirales_full_align.fasta --tree picornavirales_full_unrooted --dates picornavirales_full_treetime_date.csv --outdir picornavirales_full_unrooted_out

treetime --aln picornavirales_full_align.fasta --tree picornavirales_full_unrooted --dates picornavirales_full_treetime_date.csv --outdir picornavirales_full_unrooted_skyline_out --coalescent skyline

```

picornavirales partial
```
cd Desktop/mada-bat-picornavirus/TreeTime/picornavirales_partial

treetime --aln picornavirales_partial_align.fasta --tree picornavirales_partial_rooted --dates picornavirales_partial_treetime_date.csv --outdir picornavirales_partial_rooted_out

treetime --aln picornavirales_partial_align.fasta --tree picornavirales_partial_rooted --dates picornavirales_partial_treetime_date.csv --outdir picornavirales_partial_rooted_skyline_out --coalescent skyline

treetime --aln picornavirales_partial_align.fasta --tree picornavirales_partial_unrooted --dates picornavirales_partial_treetime_date.csv --outdir picornavirales_partial_unrooted_out

treetime --aln picornavirales_partial_align.fasta --tree picornavirales_partial_unrooted --dates picornavirales_partial_treetime_date.csv --outdir picornavirales_partial_unrooted_skyline_out --coalescent skyline

```


---
4. The R code to edit each tree is within each respective folder

Note: I think it's just because I'm using a Mac but you NEED to have the date format correctly with the csv file open when you run these commands. If you don't, it just goes back to the mm/dd/yy format which will give you incorrect values



**From the picornavirales all genome and the picornaviridae full genome tree, I think I should make 4 refseq trees for each represented family: picornaviridae, caliciviridae, secoviridae, and iflaviridae. I followed the same methods as above to do this, and the commands are as follows: 
