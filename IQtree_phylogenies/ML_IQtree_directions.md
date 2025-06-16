Building a phylogenetic tree in IQtree
---

This tutorial outlines methods for building a maximum-likelihood phylogenetic tree using the program IQtree.

---
Full genome trees for picornaviridae and caliciviridae: 
___
1. I first went to NCBI virus [here](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/) and selected the following results for each corresponding tree (I downloaded full genomes only) : 
- Picornaviridae taxid 12058 - RefSeq sequences only
- Unclassified picornaviridae taxid 478825 - RefSeq sequences only
- Caliciviridae taxid 11974 - RefSeq sequences only
- Unclassified caliciviridae taxid 179239 - RefSeq sequences only
- Kobuvirus taxid 194960
- Unclassified kobuvirus taxid 655314
- Sapelovirus taxid 686982
- Unclassified sapelovirus taxid 1073966
- Teschovirus taxid 118139
- Unclassified teschovirus taxid 2004714
- Mischivirus taxid 1511778
- Kunsagivirus taxid 1755589
- Sapovirus taxid 95341
- Unclassified sapovirus taxid 371833
- Cardiovirus taxid 12103
- Unclassified cardiovirus taxid 434308
- Hepatovirus taxid 12091
- Unclassified hepatovirus taxid 1714618
- Shanbavirus taxid 2169640
- Unclassified shanbavirus taxid 2788549
- Bat picornavirus taxid 1281456


For each genus, I downloaded a .csv file which I processed in the attached R scripts [here]() to meet the following criteria: 

-I  downloaded full nucleotide sequences (which may include RefSeq sequences), in addition to all bat sequences over 2000bp. 

I then removed any duplicates and any records which were entered twice with two accession numbers after having been reviewed to reference status. The unique records were downloaded into a  folder on my home computer in bulk from NCBI using the following code loaded into my web browser (produced from the Rscripts linked [here]()):

For all trees I used an outgroup of Sindbis virus accession NC_001547.1 since it's also short in length and not within Picornavirales
---

2. Next I added our Madagascar bat sequences using the following accession numbers corresponding to each genome.

Full genomes: 
OQ818322
PP766456
OQ818316
OQ818321
OQ818329
OQ818320
OQ818328
PP766469
PP766459
OQ818317
OQ818323
OQ818324
OQ818318
RR034B_242_fecRa_S67


Partial genomes:
PP766457
PP766455
OQ818325
OQ818319
OQ818337
PP766467
OQ818340
PP766470
PP766451
OQ818348
PP766471
PP766472
OQ818347
OQ818344
PP766460 - over 2000bp but still too short for phylogeny
PP766466 - under 2000bp
PP766473 - under 2000bp
PP766454 - under 2000bp
PP766452 - under 2000bp
PP766462 - under 2000bp
PP766461 - under 2000bp
PP766458 - under 2000bp
PP766450 - under 2000bp
PP766449 - under 2000bp
PP766474 - under 2000bp
PP766463 - under 2000bp
PP766468 - under 2000bp
OQ818342
PP766464 - under 2000bp
PP766475 - under 2000bp
PP766476 - under 2000bp
OQ818343
PP766465 - under 2000bp
OQ818345 - under 2000bp
PP766453 - under 2000bp
OQ818346 - under 2000bp
PP766477 - under 2000bp
RR034B_239_fecRa_S61
RR034B_281_fecRa_S29

While I am mainly focused on the full genome sequences, I'll probably include some of the longer partial sequences as well in order to help fill out the tree...plus some genera I only have partials for like hepatovirus and cardiovirus. Any novel sequence under 2000bp was excluded from phylogenetic analysis. 

---
I uploaded the concatenated fastas into the [MAFFT program online](https://mafft.cbrc.jp/alignment/server/) for alignment, After the alignment returns, I downloaded it as a .fasta file and saved with the same name as before, but with "_align" at the end. 

---

3. Then, I used IQtree to make the trees, which automatically performs 1000 ultrafast bootstraps. Below is an example of the command line used to do that: 

```
bin/iqtree -s .fasta -bb 1000

```
The best models for the phylogenetic trees are losted below:

```
Best model for bat picornavirus all genomes: GTR+F+R5
Best model for cardiovirus all genomes: TIM2+F+R10
Best model for hepatovirus all genomes: GTR+F+I+G4
Best model for kobuvirus all genomes: GTR+F+R10
Best model for kunsagivirus all genomes: TIM2e+R2
Best model for mischivirus all genomes: GTR+F+R4
Best model for sapelovirus all genomes: TIM2+F+R10
Best model for sapovirus all genomes: GTR+F+R10
Best model for teschovirus all genomes: TIM2+F+R9
Best model for polymerase summary tree: 
```
---

4. Once IQtree finished, I imported the resulting tree into R and built a phylogenetic tree. You can easily look at the tree in [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) before cleaning it up. The scripts for making the trees are within each respective folder linked [here](). 
