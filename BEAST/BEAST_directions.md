Building a maximum-likelihood Bayesian timetree using BEAST2
---
Resources to using BEAST2 and BEAUTi: there are many helpful tutorials online [here](https://taming-the-beast.org/tutorials/)

I am basing this workflow on previous analysis done in our lab describing novel bat Nobecoviruses [here](https://github.com/brooklabteam/Mada-Bat-CoV/blob/main/Fig4/beast-tree-instructions.md) and a henipavirus [here](https://github.com/brooklabteam/angavokely-virus/tree/c135e21f706d4c40023a5b39628c3e2b803e945c/Fig4).

---
1. I decided to make a Bayesian tree for refseq picornaviridae with all full genome picornaviridae novel malagasy picornaviruses. The alignment file and modeltest results for this alignment already exist from building the maximum likelihood tree in RAxML. However, this alignment file has an outgroup of a bat coronavirus that I will remove, as a BEAST tree should not be rooted. This alignment includes 11 full genome novel Malagasy picornaviridae sequences, one nearly full novel Malagasy hepatovirus, and 166 full genome RefSeq picornaviridae from NCBI virus taxid: 12058 and unclassified picornaviridae taxid: 478825. 

see script "pre-beast-name.R"), including the accession number and the collection date. For any sequences for which only a collection year was reported, we set the middle of the year (July 15) as the corresponding date. We then aligned these sequences using the MAFFT algorithm in Geneious, and evaluated the optimal nucleotide substitution model using ModelTest-NG.


From those results, the best model was: 









```
cd Desktop/developer/mada-bat-picornavirus/TreeTime/picornavirales

treetime --aln picornavirales_all_align.fasta --tree picornavirales_all_rooted --dates picornavirales_all_treetime_date.csv --outdir picornavirales_all_rooted_out

treetime --aln picornavirales_all_align.fasta --tree picornavirales_all_rooted --dates picornavirales_all_treetime_date.csv --outdir picornavirales_all_rooted_skyline_out --coalescent skyline

treetime --aln picornavirales_all_align.fasta --tree picornavirales_all_unrooted --dates picornavirales_all_treetime_date.csv --outdir picornavirales_all_unrooted_out

treetime --aln picornavirales_all_align.fasta --tree picornavirales_all_unrooted --dates picornavirales_all_treetime_date.csv --outdir picornavirales_all_unrooted_skyline_out --coalescent skyline

```
