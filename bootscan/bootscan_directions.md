Bootscan
---
To look at potential recombination events between novel Madagascar sequences and other described reference sequences in a plot, follow the documentation linked [here](https://sray.med.som.jhmi.edu/SCRoftware/SimPlot/online_help/hs1050.htm)
---
1. This software can only run on a PC computer. Follow the directions [here](https://sray.med.som.jhmi.edu/SCRoftware/SimPlot/) to download. 
---figures
2. Once installed, I will use the nucleotide alignments created for PySimPlot to feed into Bootscan. The analysis is very similar, it just measures % permutated trees rather than aa/nt similarity. The first sequence listed in the fasta is the reference sequence, so make sure this is always a Madagascar novel sequence. 

Based on the quality of the alignments, and based on the length of the sequences (partial sequences will likely not be analyzed), the following alignments were able to provide Bootscan output. Save the files as .csv:
- Bat picornaviruses ICTV full
- Cheravirus ICTV partial
- Hepatovirus ICTV partial
- Kobuvirus ICTV full
- Kunsagivirus ICTV full
- Mischivirus ICTV full
- Sapelovirus ICTV full
- Sapovirus ICTV and some african sequences full
- Teschovirus ICTV full


Alignment files from Geneious are in their respective labeled files within the master bootscan folder, and metadata/summary files are [here]()

---
3. After finishing populating the csv files, I made files with their data that could be read into R to make the figures (the PySimPlot R scripts can be modified to take these inputs, since the data set up and plotting scripts will be similar). The gene mapping scripts can stay the same as well.  

Bootscan plots were set to a window size of 200bp with the step size set to 20bp

In making a figure, I decided to only mark a point as a breakpoint if it was >30% permuted trees. Other papers have more stringent cutoffs, but I have seen some papers describing novel picornaviruses with a much lower cutoff, the highest breakpoint was ~57%. For this reason, the final figure for recombination only includes sapovirus, kunsagivirus, and hepatovirus.

---

The R plot to generate figures and process the data is [here]().
