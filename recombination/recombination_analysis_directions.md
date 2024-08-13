RDP4 analysis
---
To look at potential recombination events between novel Madagascar sequences and other described reference sequences in a plot, follow the documentation linked [here](http://web.cbio.uct.ac.za/~darren/RDP4Manual.pdf)
---
1. This software can only run on a PC computer. Follow the directions [here](http://web.cbio.uct.ac.za/~darren/rdp.html) to download. 

2. Once installed, I will use the nucleotide alignments created for PySimPlot to feed into RDP4, the summary sheet of alignments and sequences included per alignment is [here](https://github.com/brooklabteam/mada-bat-picornavirus/blob/main/recombination/recombination_references.xlsx).

Alignment files from Geneious are in their respective labeled files within the master recombination folder, and metadata/summary files are [here](https://github.com/brooklabteam/mada-bat-picornavirus/tree/main/recombination)

---

3. Load each alignment into RDP4 then run analysis using default settings, I saved individual .csv files of the statistics of each alignment and further summarized them into one master file within the recombination folder. Then, I saved the bootscan figues as .csv files to plot as long as bootscan was significant for recombination. 

---
4. After finishing populating the csv files, I made files with their data that could be read into R to make the figures and combines the gene map script, like what is done in the PySimPlot figure

Bootscan plots were set to a window size of 200bp with the step size set to 20bp

In making a figure, I decided to only consider the bootstrap cutoff if it was >70% as suggested by the RDP4 manual.

---

The R plot to generate figures and process the data is [here](https://github.com/brooklabteam/mada-bat-picornavirus/tree/main/recombination).
