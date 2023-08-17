Building a maximum-likelihood Bayesian timetree using BEAST2
---
Resources to using BEAST2 and BEAUTi: there are many helpful tutorials online [here](https://taming-the-beast.org/tutorials/)

I am basing this workflow on previous analysis done in our lab describing novel bat Nobecoviruses [here](https://github.com/brooklabteam/Mada-Bat-CoV/blob/main/Fig4/beast-tree-instructions.md) and a henipavirus [here](https://github.com/brooklabteam/angavokely-virus/tree/c135e21f706d4c40023a5b39628c3e2b803e945c/Fig4).

---
1. I decided to make a Bayesian tree for refseq picornaviridae with all full genome picornaviridae novel malagasy picornaviruses. The alignment file and modeltest results for this alignment already exist from building the maximum likelihood tree in RAxML. However, this alignment file has an outgroup of a bat coronavirus that I will remove, as a BEAST tree should not be rooted. This alignment includes 11 full genome novel Malagasy picornaviridae sequences, one nearly full novel Malagasy hepatovirus, and 166 full genome RefSeq picornaviridae from NCBI virus taxid: 12058 and unclassified picornaviridae taxid: 478825. 

see script "pre-beast-name.R"), including the accession number and the collection date. For any sequences for which only a collection year was reported, we set the middle of the year (July 15) as the corresponding date. We then aligned these sequences using the MAFFT algorithm in Geneious, and evaluated the optimal nucleotide substitution model using ModelTest-NG.


From those results, the best model was: 

2. Now that we have our aligned files, our metadata, and our substitution model, we can generate a .xml file in BEAUTi. Resources on picking models are [here](https://justinbagley.rbind.io/2016/10/11/setting-dna-substitution-models-beast/),[here](http://www.iqtree.org/doc/Substitution-Models). 

To prepare the .xml file, we used the following parameters in tab inputs at the top of the screen in BEAUti:

- Tip Dates: We used the date of sample collection as the "Tip Date." For any sample from GenBank which only listed year of collection, we set the tip date to July 15 of the year of collection. Each alignment was uploaded to BEAUti with sequence names arranged so as to easily decipher the date.
- Site Model: Following output from ModelTest-NG, we selected a "Gamma Site Model" with Gamma Category 4 and estimated the proportion of invariate sites.
- Clock Model: We built comparative phylogenies using both a strict molecular clock and a relaxed exponential molecular clock, specified within a non-parameteric Bayesian Skyline Coalescent model, following previous approaches for bat coronvirus analyses (i.e. Lau et al. 2020.
- Priors: We used a Bayesian Skyline Coalescent model. The clock rate prior was set to a lognormal distribution with a mean of 0.001, following published values for RNA viruses (Jenkins et al. 2014), and all other priors were left at default values specified in BEAUti. Both xml files for the strick and relaxed molecular clocks are available in the subfolder "3-BEAST-Nobeco" within the Fig4 folder.
- MCMC: We used an MCMC chain length of 700,000,000 iterations and set tracelog and treelog every 10,000 iterations. All other MCMC metrics were left at default.
- Population Size: The Bayesian Skyline Coalescent model by default assumes that the population size of the dataset will change 5 times spaced evenly across the course of this sampling period. Because our available samples were limited and spanned a wide geographic area, we edited this parameter to permit only one population size. You can make this edit in BEAUti by clicking "View" in the top panel, selecting "Show Initialization panel" and then specifying the dimension of "bPopSizes" and "GroupSizes" both to be 1 instead of 5.




```
cd Desktop/developer/mada-bat-picornavirus/TreeTime/picornavirales

treetime --aln picornavirales_all_align.fasta --tree picornavirales_all_rooted --dates picornavirales_all_treetime_date.csv --outdir picornavirales_all_rooted_out

treetime --aln picornavirales_all_align.fasta --tree picornavirales_all_rooted --dates picornavirales_all_treetime_date.csv --outdir picornavirales_all_rooted_skyline_out --coalescent skyline

treetime --aln picornavirales_all_align.fasta --tree picornavirales_all_unrooted --dates picornavirales_all_treetime_date.csv --outdir picornavirales_all_unrooted_out

treetime --aln picornavirales_all_align.fasta --tree picornavirales_all_unrooted --dates picornavirales_all_treetime_date.csv --outdir picornavirales_all_unrooted_skyline_out --coalescent skyline

```
