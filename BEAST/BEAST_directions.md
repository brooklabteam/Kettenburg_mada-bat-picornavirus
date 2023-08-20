Building a maximum-likelihood Bayesian timetree using BEAST2
---
Resources to using BEAST2 and BEAUTi: there are many helpful tutorials online [here](https://taming-the-beast.org/tutorials/)

I am basing this workflow on previous analysis done in our lab describing novel bat Nobecoviruses [here](https://github.com/brooklabteam/Mada-Bat-CoV/blob/main/Fig4/beast-tree-instructions.md) and a henipavirus [here](https://github.com/brooklabteam/angavokely-virus/tree/c135e21f706d4c40023a5b39628c3e2b803e945c/Fig4).

---
1. I decided to make a Bayesian tree for refseq picornaviridae with all full genome picornaviridae novel malagasy picornaviruses. The alignment file and modeltest results for this alignment already exist from building the maximum likelihood tree in RAxML. However, this alignment file has an outgroup of a bat coronavirus that I will remove, as a BEAST tree should not be rooted. This alignment includes 11 full genome novel Malagasy picornaviridae sequences,and 166 full genome RefSeq picornaviridae from NCBI virus taxid: 12058 and unclassified picornaviridae taxid: 478825. 

see script "pre-beast-name.R"), including the accession number and the collection date. For any sequences for which only a collection year was reported, we set the middle of the year (July 15) as the corresponding date. We then aligned these sequences using the MAFFT algorithm in Geneious, and evaluated the optimal nucleotide substitution model using ModelTest-NG.

From those results, the best model was: GTR+I+G4

lnL:                -1319332.0078
Frequencies:        0.2552 0.2348 0.2250 0.2850
Subst. Rates:       1.7049 3.5861 1.4533 1.3349 3.6795 1.0000 
Inv. sites prop:    0.0211
Gamma shape:        1.8944
Score:              2639390.0156
Weight:             0.8383

For caliciviridae refseq tree (53 full refseq genomes) with full genome novel malagasy sapovirus (1 genome)
Best model: TVM+I+G4

lnL:                -361848.4774
Frequencies:        0.2373 0.2795 0.2427 0.2404
Subst. Rates:       2.0227 3.8404 1.6143 1.1945 3.8404 1.0000 
Inv. sites prop:    0.0220
Gamma shape:        1.4559
Score:              723921.9548
Weight:             0.8212



2. Now that we have our aligned files, our metadata, and our substitution model, we can generate a .xml file in BEAUTi. Resources on picking models are [here](https://justinbagley.rbind.io/2016/10/11/setting-dna-substitution-models-beast/),[here](http://www.iqtree.org/doc/Substitution-Models). 

To prepare the .xml file, we used the following parameters in tab inputs at the top of the screen in BEAUti:

- Tip Dates: We used the date of sample collection as the "Tip Date." For any sample from GenBank which only listed year of collection, we set the tip date to July 15 of the year of collection. Each alignment was uploaded to BEAUti with sequence names arranged so as to easily decipher the date.

- Site Model: Following output from ModelTest-NG, we selected a "Gamma Site Model" with Gamma category 4 and an estimated 0.001 for the proportion invariant. For caliciviridae, TVM is a less common model, so according to the blog post [here](https://justinbagley.rbind.io/2016/10/11/setting-dna-substitution-models-beast/) I left it as GTR but unchecked the AG rate operator. 

- Clock Model: following literature on other picornaviridae bayesian trees, I opted to use both a relaxed lognormal molecular clock and a strict molecular clock. Most papers use a relaxed clock, but some picornaviridae genera perform better under strict clocks, see [bayesian_paper_refereces]() for my summaries of different papers.

- Priors: We used a Bayesian Skyline Coalescent model for both relaxed and strict clocks. The clock rate prior was set based on the comparisons done in [Hicks et al 2011](), I set the number of substitutions to 1.60 × 10−3 ns/s/y (0.0016) which was the average of non-enterovirus rates, and all other priors were left at default values specified in BEAUti. Both xml files for the strict and relaxed molecular clocks are available in the folder [here]().
  Note: for Sapovirus the average substitution rate used in a caliciviridae refseq bayesian tree   was 0.00235, based on [Tohma et al. 2020]()

- MCMC:  We used an MCMC chain length of 100,000,000 iterations and set tracelog and treelog every 10,000 iterations, with a 10% burn-in. All other MCMC metrics were left at default.

- Population Size: The Bayesian Skyline Coalescent model by default assumes that the population size of the dataset will change 5 times spaced evenly across the course of this sampling period. Because our available samples were limited and spanned a wide geographic area, we edited this parameter to permit only one population size. You can make this edit in BEAUti by clicking "View" in the top panel, selecting "Show Initialization panel" and then specifying the dimension of "bPopSizes" and "GroupSizes" both to be 1 instead of 5.


3. Visualizing the tree: The initial 10% of MCMC iterations were removed as burn-in. Parameter convergence was assessed visually using Tracer v1.6. Once verified, we used TreeAnnotator to average across the BEAST tree output, then visualized the resulting tree in R using the script [here]().
