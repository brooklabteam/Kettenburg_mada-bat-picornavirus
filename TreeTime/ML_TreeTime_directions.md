Building a maximum-likelihood TimeTree using RAxML output
---

This tutorial outlines methods for building a maximum-likelihood time-scaled phylogenies described [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5758920/) with tutorials listed [here](https://github.com/neherlab/treetime). I will use the same tree files as done in building the ML trees. Each is in their own folder in the master folder [here](). I chose this method over BEAST as my picornavirales dataset would take an extremely long time to run. 

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
cd Desktop/developer/mada-bat-picornavirus/TreeTime/picornavirales

treetime --aln picornavirales_all_align.fasta --tree picornavirales_all_rooted --dates picornavirales_all_treetime_date.csv --outdir picornavirales_all_rooted_out

treetime --aln picornavirales_all_align.fasta --tree picornavirales_all_rooted --dates picornavirales_all_treetime_date.csv --outdir picornavirales_all_rooted_skyline_out --coalescent skyline

treetime --aln picornavirales_all_align.fasta --tree picornavirales_all_unrooted --dates picornavirales_all_treetime_date.csv --outdir picornavirales_all_unrooted_out

treetime --aln picornavirales_all_align.fasta --tree picornavirales_all_unrooted --dates picornavirales_all_treetime_date.csv --outdir picornavirales_all_unrooted_skyline_out --coalescent skyline

```

bat picornavirales
```
cd Desktop/developer/mada-bat-picornavirus/TreeTime/bat_picornavirales

treetime --aln bat_picornavirales_all_over_3kb_align.fasta --tree bat_picornavirales_all_rooted --dates bat_picornavirales_all_dates.csv --outdir bat_picornavirales_all_rooted_out

treetime --aln bat_picornavirales_all_over_3kb_align.fasta --tree bat_picornavirales_all_rooted --dates bat_picornavirales_all_dates.csv --outdir bat_picornavirales_all_rooted_skyline_out --coalescent skyline

treetime --aln bat_picornavirales_all_over_3kb_align.fasta --tree bat_picornavirales_all_unrooted --dates bat_picornavirales_all_dates.csv --outdir bat_picornavirales_all_unrooted_out

treetime --aln bat_picornavirales_all_over_3kb_align.fasta --tree bat_picornavirales_all_unrooted --dates bat_picornavirales_all_dates.csv --outdir bat_picornavirales_all_unrooted_skyline_out --coalescent skyline

```

caliciviridae
```
cd Desktop/developer/mada-bat-picornavirus/TreeTime/caliciviridae

treetime --aln caliciviridae_refseq_align.fasta --tree caliciviridae_rooted --dates caliciviridae_dates.csv --outdir caliciviridae_rooted_out

treetime --aln caliciviridae_refseq_align.fasta --tree caliciviridae_rooted --dates caliciviridae_dates.csv --outdir caliciviridae_rooted_skyline_out --coalescent skyline

treetime --aln caliciviridae_refseq_align.fasta --tree caliciviridae_unrooted --dates caliciviridae_dates.csv --outdir caliciviridae_all_unrooted_out

treetime --aln caliciviridae_refseq_align.fasta --tree caliciviridae_unrooted --dates caliciviridae_dates.csv --outdir caliciviridae_unrooted_skyline_out --coalescent skyline
```

iflaviridae
```
cd Desktop/developer/mada-bat-picornavirus/TreeTime/iflaviridae

treetime --aln iflaviridae_refseq_align.fasta --tree iflaviridae_rooted --dates iflaviridae_dates.csv --outdir iflaviridae_rooted_out

treetime --aln iflaviridae_refseq_align.fasta --tree iflaviridae_rooted --dates iflaviridae_dates.csv --outdir iflaviridae_rooted_skyline_out --coalescent skyline

treetime --aln iflaviridae_refseq_align.fasta --tree iflaviridae_unrooted --dates iflaviridae_dates.csv --outdir iflaviridae_all_unrooted_out

treetime --aln iflaviridae_refseq_align.fasta --tree iflaviridae_unrooted --dates iflaviridae_dates.csv --outdir iflaviridae_unrooted_skyline_out --coalescent skyline
```

picornaviridae
```
cd Desktop/developer/mada-bat-picornavirus/TreeTime/picornaviridae

treetime --aln picornaviridae_refseq_align.fasta --tree picornaviridae_rooted --dates picornaviridae_dates.csv --outdir picornaviridae_rooted_out

treetime --aln picornaviridae_refseq_align.fasta --tree picornaviridae_rooted --dates picornaviridae_dates.csv --outdir picornaviridae_rooted_skyline_out --coalescent skyline

treetime --aln picornaviridae_refseq_align.fasta --tree picornaviridae_unrooted --dates picornaviridae_dates.csv --outdir picornaviridae_all_unrooted_out

treetime --aln picornaviridae_refseq_align.fasta --tree picornaviridae_unrooted --dates picornaviridae_dates.csv --outdir picornaviridae_unrooted_skyline_out --coalescent skyline
```

secoviridae
```
cd Desktop/developer/mada-bat-picornavirus/TreeTime/secoviridae

treetime --aln secoviridae_refseq_align.fasta --tree secoviridae_rooted --dates secoviridae_dates.csv --outdir secoviridae_rooted_out

treetime --aln secoviridae_refseq_align.fasta --tree secoviridae_rooted --dates secoviridae_dates.csv --outdir secoviridae_rooted_skyline_out --coalescent skyline

treetime --aln secoviridae_refseq_align.fasta --tree secoviridae_unrooted --dates secoviridae_dates.csv --outdir secoviridae_all_unrooted_out

treetime --aln secoviridae_refseq_align.fasta --tree secoviridae_unrooted --dates secoviridae_dates.csv --outdir secoviridae_unrooted_skyline_out --coalescent skyline
```


---
4. The R code to edit each tree is within each respective folder

Note: I think it's just because I'm using a Mac but you NEED to have the date format correctly with the csv file open when you run these commands. If you don't, it just goes back to the mm/dd/yy format which will give you incorrect values

