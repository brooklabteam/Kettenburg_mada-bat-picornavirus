PySimPlot
---
To look at the similarity between novel Madagascar sequences and other described reference sequences in a plot, follow the documentation linked [here](https://github.com/jonathanrd/pySimPlot)
---
1. To install on your local computer (MacOS for me), open the command line and type

```
pip3 install pysimplot
```
You may have to update your version of Python if you have a Mac M1 chip

---
2. Once installed, you will need to align your novel sequences with the sequences you wish to compare it against. The first sequence listed in the fasta will be what PySimPlot is comparing everything against so depending on how you want to present the data, you may want to list your novel sequences first or you may want to list your reference sequence(s) first. I just aligned everything in Geneious but you can use [MAFFT](https://mafft.cbrc.jp/alignment/software/) again too. 


I have an excel sheet summarizing the comparisons I picked for PySimPlot similarity analysis and pairwise BLAST comparisons [here]().

I want to make a figure of the following: 
- Pairwise comparison with each respective virus's top aa BLAST hit (noit going into figure, may use distance matrix in a figure)
- Representative sequences for species defined in ICTV and the top BLAST hit (going into figure)
- A matching virus genus if it has been found in the Africa/SWIO region in a bat...want to show how similar/different these Malagasy bat viruses are those to the picornaviruses found in African bats (either a table or a distance matrix will be used in the final figure)

Malagasy sequences included in PySimplot figure (ICTV+BLAST sequences)
OQ818316	mischivirus
OQ818317	kunsagivirus
OQ818318	teschovirus
OQ818320	sapelovirus
OQ818321	sapelovirus
OQ818322	kobuvirus
OQ818323	teschovirus
OQ818324	teschovirus
OQ818329	sapelovirus
OQ818337	hepatovirus
OQ818344	sapelovirus

Malagasy sequences used in African bat comparison table
OQ818337	hepatovirus
OQ818316	mischivirus
OQ818317	kunsagivirus
OQ818319	sapovirus
OQ818320	sapelovirus
OQ818321	sapelovirus
OQ818322	kobuvirus
OQ818329	sapelovirus
OQ818335	felisavirus
OQ818340	sapovirus
OQ818341	felisavirus
OQ818342	sapelovirus
OQ818343	sapelovirus
OQ818344	sapelovirus
OQ818345	sapovirus
OQ818347	sapovirus
OQ818348	sapovirus

All other sequences (the aa sequences not included in Fig 2 and all nt versions of the figures) will be included in a supplemental figure. There is also a supplemental figure of the african bat picornavirus pysimplots. 

Alignment and distance matrix output from Geneious are in their respective labeled files within the master PySimPlot folder, and metadata/summary files are [here]()

---
3. Following the template below, I made similarity plots:

cd into each folder where your aligned fasta files are before running this command:

```
pysimplot -i aligned.fasta -o data.csv

```

Where aligned.fasta is your aligned fasta file and data.csv is something you name for pysimplot to output data to. A blank csv file does not need to be created beforehand.

---
4. After finishing populating the csv files, I made files with their data that could be read into R to make the figures. 

Amino acid plots were set to 100aa window size with the step size of 20aa, nucleotide plots were set to a window size of 200bp with the step size set to 20bp
---

The R plot to generate figures and process the data is [here]().
