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

I have an excel sheet summarizing the comparisons I picked for PySimPlot similarity analysis. They include full-genome novel sequences, the top hits via BLASTx, and those that phylogenetically clustered with the novel sequences [here](https://github.com/brooklabteam/mada-bat-picornavirus/blob/main/PySimPlot/pysimplot_references.xlsx).

All other sequences (the aa sequences not included in Fig 2 and all nt versions of the figures) will be included in a supplemental figure. 

All aligned files are [here](https://github.com/brooklabteam/mada-bat-picornavirus/tree/main/PySimPlot/fastas)

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

The R plot to generate figures and process the data is [here](https://github.com/brooklabteam/mada-bat-picornavirus/blob/main/PySimPlot/pysimplot_plotting.R).
