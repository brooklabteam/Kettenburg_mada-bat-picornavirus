Building a phylogenetic tree in RAxML
---

This tutorial outlines methods for building a maximum-likelihood phylogenetic tree using the program RAxML.

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


For each picorna genus, I downloaded a .csv file which I processed in the attached R scripts [here]() to meet the following criteria: 
- For RefSeq full genome tree: all RefSeq genome picornas with a bat host, all RefSeq genome picorna reference sequences with a non-bat host, and all RefSeq genome unclassified picornas that grouped with the picorna genera with a bat host.

-For all other full genome trees, I just downloaded full nucleotide sequences (which may include RefSeq sequences) and repeated the same criteria.

I then removed any duplicates and any records which were entered twice with two accession numbers after having been reviewed to reference status. The unique records were downloaded into a  folder on my home computer in bulk from NCBI using the following code loaded into my web browser (produced from the Rscripts linked [here]()):

For all trees I used an outgroup of Sindbis virus accession NC_001547.1 since it's also short in length and not within Picornavirales
---

2. Next I added our Madagascar bat sequences using the following accession numbers corresponding to each genome.

Full genomes: 

Partial genomes: 

---
I uploaded the concatenated fastas into the [MAFFT program online](https://mafft.cbrc.jp/alignment/server/) for alignment, After the alignment returns, I downloaded it as a .fasta file and saved with the same name as before, but with "_align" at the end. I then used the pre-prep R file used previously to download accession numbers from NCBI [here]() to edit the names of each sequence in the MSA, since RAxML won't accept spaces, periods, dashes, slashes, colon, semicolons, or parentheses in the sequence names.

However, if you do the alignments in Geneious, you can export and select the option to export the fasta without the sequence description, so you would not have to do the above cleaning step in the R script. 

---

5. Next, I used the MSA to compare nucleotide substitution models in the program [ModelTest-NG](https://github.com/ddarriba/modeltest). This can be run on the command line on your personal computer but it will take awhile. I ran them on Midway instead, see Brook lab documentation on how to do this [here](https://github.com/brooklabteam/brooklab-resources/blob/main/modeltest-ng.md). See the same repo on how to install RAxML on the cluster [here](https://github.com/brooklabteam/brooklab-resources/blob/main/RAxML-mpi.md).

Below are the batch scripts I used for ModelTest for picornaviridae and picornavirales: 

I have trouble putting the module load and conda step within the batch script, so just run this before submitting scripts: 

```
module load python
conda activate /project2/cbrook/software/modeltest_env
```

ModelTest for all (partial and full in one tree) picornavirales
```
#!/bin/bash
#SBATCH --job-name=picornavirales_all
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks=8
#SBATCH --time=36:00:00

module load python
conda activate /project2/cbrook/software/modeltest_env

modeltest-ng -i picornavirales_all_align.fasta -d nt -t ml -p 8
```

ModelTest for full picornavirales
```
#!/bin/bash
#SBATCH --job-name=picornavirales_full
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks=8
#SBATCH --time=36:00:00

module load python
conda activate /project2/cbrook/software/modeltest_env

modeltest-ng -i picornavirales_full_align.fasta -d nt -t ml -p 8
```

ModelTest for partial picornavirales
```
#!/bin/bash
#SBATCH --job-name=picornavirales_partial
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks=8
#SBATCH --time=36:00:00

module load python
conda activate /project2/cbrook/software/modeltest_env

modeltest-ng -i picornavirales_partial_align.fasta -d nt -t ml -p 8
```

ModelTest for all (partial and full in one tree) picornaviridae
```
#!/bin/bash
#SBATCH --job-name=picornaviridae_all
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks=8
#SBATCH --time=36:00:00

module load python
conda activate /project2/cbrook/software/modeltest_env

modeltest-ng -i picornaviridae_all_align.fasta -d nt -t ml -p 8
```

ModelTest for full picornaviridae
```
#!/bin/bash
#SBATCH --job-name=picornaviridae_full
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks=8
#SBATCH --time=36:00:00

module load python
conda activate /project2/cbrook/software/modeltest_env

modeltest-ng -i picornaviridae_full_align.fasta -d nt -t ml -p 8
```

ModelTest for partial picornaviridae
```
#!/bin/bash
#SBATCH --job-name=picornaviridae_partial
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks=8
#SBATCH --time=36:00:00

module load python
conda activate /project2/cbrook/software/modeltest_env

modeltest-ng -i picornaviridae_partial_align.fasta -d nt -t ml -p 8
```

I saved the results of ModelTest-NG [here](). 

Best model for picornavirales full genomes: GTR+G4
Best model for picornavirales partial genomes: GTR+G4
Best model for picornavirales all genomes: GTR+G4

Best model for picornaviridae full genomes: GTR+G4
Best model for picornaviridae partial genomes: GTR+G4
Best model for picornaviridae all genomes: TVM+G4

---

6. Once ModelTest-NG finishes (it took about 6 days for the picornavirales and maybe 2 days for picornaviridae), it is time to build a maximum likelihood tree using RAxML. See documentation on their website for how to get this running on your home computer and/or computing cluster. 

Check whether the alignment can be read by RAxML and the number of threads you should use (I did these on my local computer so the command is raxml-ng, if you run these on the cluster be sure to use raxml-ng-mpi:

Picornavirales all (partial and full) genomes
```
raxml-ng --check --msa picornavirales_all_align.fasta --model GTR+G4 --prefix T100
raxml-ng --parse --msa picornavirales_all_align.fasta --model GTR+G4 --prefix T101
```
Use 12 threads

Picornavirales full genomes
```
raxml-ng --check --msa picornavirales_full_align.fasta --model GTR+G4 --prefix T102
raxml-ng --parse --msa picornavirales_full_align.fasta --model GTR+G4 --prefix T103
```
Use 12 threads

Picornavirales partial genomes
```
raxml-ng --check --msa picornavirales_partial_align.fasta --model GTR+G4 --prefix T104
raxml-ng --parse --msa picornavirales_partial_align.fasta --model GTR+G4 --prefix T105
```
Use 12 threads

Picornaviridae all (partial and full) genomes
```
raxml-ng --check --msa picornaviridae_all_align.fasta --model TVM+G4 --prefix T106
raxml-ng --parse --msa picornaviridae_all_align.fasta --model TVM+G4 --prefix T107
```
Use 10 threads

Picornaviridae full genomes
```
raxml-ng --check --msa picornaviridae_full_align.fasta --model GTR+G4 --prefix T108
raxml-ng --parse --msa picornaviridae_full_align.fasta --model GTR+G4 --prefix T109
```
Use 9 threads

Picornaviridae partial genomes
```
raxml-ng --check --msa picornaviridae_partial_align.fasta --model GTR+G4 --prefix T110
raxml-ng --parse --msa picornaviridae_partial_align.fasta --model GTR+G4 --prefix T111
```
Use 10 threads

---

7. Finally, I kicked off RAxML in each subfolder using the SLURM scripts using below, this does the ML tree search and bootstrapping, but you can do the processes separately to try and save time if the cluster times out on the run time, see [here](https://github.com/amkozlov/raxml-ng/wiki/Tutorial) for more resources on running RAxML:

The script below are for a one of command that does the tree inference, bootstrapping, and calculated support. If the picornavirales trees time out, I'll do the steps separately. 

Picornavirales all (partial and full) genomes
```
#!/bin/bash
#SBATCH --job-name=picornavirales_all
#SBATCH --partition=broadwl
#SBATCH --output=picornavirales_all.out
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --ntasks-per-node=12
#SBATCH --time=36:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1

raxml-ng-mpi --all --msa picornavirales_all_align.fasta --model GTR+G4 --prefix T1  --seed 1 --threads auto{12} --bs-metric fbp,tbe
```


Picornavirales full genomes
```
#!/bin/bash
#SBATCH --job-name=picornavirales_full
#SBATCH --partition=broadwl
#SBATCH --output=picornavirales_full.out
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --ntasks-per-node=12
#SBATCH --time=36:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1

raxml-ng-mpi --all --msa picornavirales_full_align.fasta --model GTR+G4 --prefix T1  --seed 1 --threads auto{12} --bs-metric fbp,tbe
```


Picornavirales partial genomes
```
#!/bin/bash
#SBATCH --job-name=picornavirales_partial
#SBATCH --partition=broadwl
#SBATCH --output=picornavirales_partial.out
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --ntasks-per-node=12
#SBATCH --time=36:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1

raxml-ng-mpi --all --msa picornavirales_partial_align.fasta --model GTR+G4 --prefix T1  --seed 1 --threads auto{12} --bs-metric fbp,tbe
```


Picornaviridae all (partial and full) genomes
```
#!/bin/bash
#SBATCH --job-name=picornaviridae_all
#SBATCH --partition=broadwl
#SBATCH --output=picornaviridae_all.out
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --ntasks-per-node=10
#SBATCH --time=36:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1

raxml-ng-mpi --all --msa picornaviridae_all_align.fasta --model TVM+G4 --prefix T1  --seed 1 --threads auto{10} --bs-metric fbp,tbe
```


Picornaviridae full genomes
```
#!/bin/bash
#SBATCH --job-name=picornaviridae_full
#SBATCH --partition=broadwl
#SBATCH --output=picornaviridae_full.out
#SBATCH --nodes=1
#SBATCH --ntasks=9
#SBATCH --ntasks-per-node=9
#SBATCH --time=36:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1

raxml-ng-mpi --all --msa picornaviridae_full_align.fasta --model GTR+G4 --prefix T1  --seed 1 --threads auto{9} --bs-metric fbp,tbe
```


Picornaviridae partial genomes
```
#!/bin/bash
#SBATCH --job-name=picornaviridae_partial
#SBATCH --partition=broadwl
#SBATCH --output=picornaviridae_partial.out
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --ntasks-per-node=10
#SBATCH --time=36:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1

raxml-ng-mpi --all --msa picornaviridae_partial_align.fasta --model GTR+G4 --prefix T1  --seed 1 --threads auto{10} --bs-metric fbp,tbe
```

---

8. Once RAxML finished, I imported the resulting tree into R and built a phylogenetic tree. You can easily look at the tree in [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) before cleaning it up. The scripts for making the trees are within each respective folder. 

---

9. Upon viewing the picornavirales and the picornaviridae trees, I decided to make additonal trees for each virus family that I believe the novel Madagascar sequences belong to. So I'll have a picornavirales, a picornaviridae, a secoviridae, a caliciviridae, and an iflaviridae tree. To get the refseq's for these, I followed the same instructions as above using the same outgroup. 

Also following the same instructions as above and using the same outgroup, I decided to make a non-refseq picornavirales tree with only sequences from a bat host. Once I downloaded the sequences from NCBI virus, I put them into Geneious to check for duplicate sequences (they may be different viruses but if the sequences are super similar RAxmL and modeltest give me a hard time), and then removed them from the file before aligning. With the bat picornavirales, I tested a variety of different alignments and ended up choosing the following to use for bat picornavirales full reference sequences and bat picornavirales all reference sequences
- bat picornavirales full reference sequences with mada sequences >3kb, as done in the picornavirales refseq tree
- bat picornavirales all reference sequences with mada and reference sequences both over 3kb

The commands I used in modeltest and RAxML are listed below. 
---

ModelTest bash scripts: 

Bat picornavirales full reference sequences all over 3kb
```
#!/bin/bash
#SBATCH --job-name=bat_picorna_all
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks=8
#SBATCH --time=36:00:00

module load python
conda activate /project2/cbrook/software/modeltest_env

modeltest-ng -i bat_picornavirales_all_over_3kb_align.fasta -d nt -t ml -p 8
```
Best model: GTR+G4 
Number of threads for RAxmL:8



Bat picornavirales all reference sequences all over 3kb
```
#!/bin/bash
#SBATCH --job-name=batpicorna_full
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks=8
#SBATCH --time=36:00:00

module load python
conda activate /project2/cbrook/software/modeltest_env

modeltest-ng -i bat_picornavirales_full_3kb_align.fasta -d nt -t ml -p 8
```
Best model: GTR+G4
Number of threads for RAxmL:8



Caliciviridae refseq
```
#!/bin/bash
#SBATCH --job-name=cal
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks=8
#SBATCH --time=36:00:00

module load python
conda activate /project2/cbrook/software/modeltest_env

modeltest-ng -i caliciviridae_refseq_align.fasta -d nt -t ml -p 8
```
Best model: TVM+I+G4
Number of threads for RAxmL: 8



Iflaviridae refseq
```
#!/bin/bash
#SBATCH --job-name=ifla
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks=8
#SBATCH --time=36:00:00

module load python
conda activate /project2/cbrook/software/modeltest_env

modeltest-ng -i iflaviridae_refseq_align.fasta -d nt -t ml -p 8
```
Best model: GTR+G4
Number of threads for RAxmL: 9



Picornaviridae refseq
```
#!/bin/bash
#SBATCH --job-name=picornaviridae_all
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks=8
#SBATCH --time=36:00:00

module load python
conda activate /project2/cbrook/software/modeltest_env

modeltest-ng -i picornaviridae_refseq_align.fasta -d nt -t ml -p 8
```
Best model: GTR+G4
Number of threads for RAxmL: 9



Secoviridae refseq
```
#!/bin/bash
#SBATCH --job-name=seco
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks=8
#SBATCH --time=36:00:00

module load python
conda activate /project2/cbrook/software/modeltest_env

modeltest-ng -i secoviridae_refseq_align.fasta -d nt -t ml -p 8
```
Best model: GTR+G4
Number of threads for RAxmL: 8





---
Scripts for RAxML trees

Bat picornavirales all over 3kb
```
#!/bin/bash
#SBATCH --job-name=batpicorna_all
#SBATCH --partition=broadwl
#SBATCH --output=batpicorna_all.out
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --time=36:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1

raxml-ng-mpi --all --msa bat_picornavirales_all_over_3kb_align.fasta --model GTR+G4 --prefix T1  --seed 1 --threads auto{8} --bs-metric fbp,tbe
```


Bat picornavirales full
```
#!/bin/bash
#SBATCH --job-name=picornavirales_all
#SBATCH --partition=broadwl
#SBATCH --output=picornavirales_all.out
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --time=36:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1

raxml-ng-mpi --all --msa bat_picornavirales_full_3kb_align.fasta --model GTR+G4 --prefix T1  --seed 1 --threads auto{8} --bs-metric fbp,tbe
```


Secoviridae refseq
```
#!/bin/bash
#SBATCH --job-name=seco
#SBATCH --partition=broadwl
#SBATCH --output=seco.out
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --time=36:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1

raxml-ng-mpi --all --msa secoviridae_refseq_align.fasta --model GTR+G4 --prefix T1  --seed 1 --threads auto{8} --bs-metric fbp,tbe
```


Picornaviridae refseq
```
#!/bin/bash
#SBATCH --job-name=picornaviridae
#SBATCH --partition=broadwl
#SBATCH --output=picornaviridae.out
#SBATCH --nodes=1
#SBATCH --ntasks=9
#SBATCH --ntasks-per-node=9
#SBATCH --time=36:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1

raxml-ng-mpi --all --msa picornaviridae_refseq_align.fasta --model GTR+G4 --prefix T1  --seed 1 --threads auto{9} --bs-metric fbp,tbe

```


Iflaviridae refseq
```
#!/bin/bash
#SBATCH --job-name=ifla
#SBATCH --partition=broadwl
#SBATCH --output=ifla.out
#SBATCH --nodes=1
#SBATCH --ntasks=9
#SBATCH --ntasks-per-node=9
#SBATCH --time=36:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1

raxml-ng-mpi --all --msa iflaviridae_refseq_align.fasta --model GTR+G4 --prefix T1  --seed 1 --threads auto{9} --bs-metric fbp,tbe
```


Caliciviridae refseq
```
#!/bin/bash
#SBATCH --job-name=cal
#SBATCH --partition=broadwl
#SBATCH --output=call.out
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --time=36:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1

raxml-ng-mpi --all --msa caliciviridae_refseq_align.fasta --model TVM+I+G4 --prefix T1  --seed 1 --threads auto{8} --bs-metric fbp,tbe
```
---

