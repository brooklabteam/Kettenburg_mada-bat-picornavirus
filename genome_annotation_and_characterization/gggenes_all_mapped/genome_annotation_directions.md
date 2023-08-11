Using gggenes to make genome annotation figures
---

This tutorial outlines methods for making figures in R of full and partial genomes using the R package gggenes. I outline directions below.

---

For gggenes figure: 

First, I downloaded the sequences in Geneious using the NCBI accession numbers to make sure all the annotations were up to date. From there, I exported a csv files of all the annotations of the full sequences and the partial sequences. From there, I made csv files for all full sequences, all partial sequences, all full and partial by each virus type, and with and without polyprotein so I can show the location of all the peptides. I made each of these their own csv file in case I needed them, but you can subset the data in R so you don't need separate files. 

In these files you need the following information: 
 - A column named "molecule" which has the accession numbers
 - A column named "start" which has the starting genome position for the gene
 - A column named "end" which has the ending genome position for the gene
 - A column named "gene" which as each annotation name
 - A column named "virus" only for the files that have all the sequences together so you can keep track of them
 
 
Lastly, view my R scripts to see the steps I followed to make the figured but you can also view the GitHub pages for gggenes [here] 