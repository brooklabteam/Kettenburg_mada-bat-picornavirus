Mapping reads for full and partial sequences
---

This tutorial outlines methods for mapping the reads for contigs pulled from CZID in order to check coverage. I outline directions below.

---

At this point, I have a folder in Geneious with all of the hiqual sequences I pulled out of the command line. I made a duplicate of it and labeled it as read mapping so it's all contained in one folder.  

Next you should go to CZID [here]() and go to the respective folder for the sample you're looking for, I have most of my samples from the RR034B1_feces folder, but a couple from RR034B2_urine_wholeblood_novaseq as well. Once you go into each folder, select a specific sample that needs a read mapped. 

After clicking on a sample, go to the upper right corner of the screen to the "Download button", click "that "View results folder". Scroll all the way to the bottom to "Step Nonhost Fastq: Filter original fastq/fasta input files to only contain non-host reads processed by CZ ID". Download both files. They should have .R1 and .R2 at the end of them. 

Once these are downloaded, drag the files to the read mapping folder in Geneious. Geneious will ask you to confirm the sequencing technology and will ask whether you want to pair the reads. Pair them and select Illumina as the sequencing technology. 

Now you have your non-host reads for a given sample. Select the non-host reads file and your sequence of interest and click align/assemble -> map to reference. Make sure your sequence is the sequence of interest for the reads to map to. I left all the default settings as is then mapped the reads. 


In the resulting file, it will show you the coverage per position in the genome. We're looking for decent coverage of at least two reads per site, ideally with overlapping positions to cover all aread of the genome. If the sequence does not have sufficient coverage, I make a file labeled "discarded" and put the sequence and the corresponding read report files in there. These sequences will not be further analyzed. On those that do have good coverage, click the graphs tab on the bottom right and select "export" to get a csv file of the coverage per position. From this, you can make one per sequence and use these to make graphs in R. 


These values are "coverage" and we want the data in reads per million. We can get the total reads from CZID, per sample, go to the upper right corner and click "sample info". Then, click the tab that says pipelines and you'll see the total reads value. That reads value should be divided by the scaling factor of (total reads overall/1,000,000). The resulting value will be used in the csv sheet from Geneious. Divide every genome position value by this scaling factor to transform coverage to reads per million. 


Full genome sequences (accession number and CZID ID): 
-OQ818316/RR034B_029_NODE_84     62,438,032 total reads/1,000,000 = 62.438 rPM
-OQ818317/RR034B_079_NODE_3      117,276,390 total reads/1,000,000 = 117.276 rPM
-OQ818318/RR034B_094_NODE_3      32,858,460 total reads/1,000,000 = 32.858 rPM
-OQ818319/RR034B_096_NODE_2      44,295,222 total reads/1,000,000 = 44.295 rPM
-OQ818320/RR034B_163_NODE_204    46,343,504 total reads/1,000,000 = 46.343 rPM
-OQ818321/RR034B_189_NODE_1      32,464,964 total reads/1,000,000 = 32.464 rPM
-OQ818322/RR034B_422_NODE_1      17,718,526 total reads/1,000,000 = 17.718 rPM
-OQ818323/RR034B_244_NODE_40     95,323,590 total reads/1,000,000 = 95.323 rPM
-OQ818324/RR034B_259_NODE_68     18,588,814 total reads/1,000,000 = 18.588 rPM
-OQ818325/RR034B_268_NODE_77     31,773,930 total reads/1,000,000 = 31.773 rPM
-OQ818326/RR034B_279_NODE_8      32,374,956 total reads/1,000,000 = 32.374 rPM
-OQ818327/RR034B_279_NODE_16     32,374,956 total reads/1,000,000 = 32.374 rPM
-OQ818328/RR034B_288_NODE_100    150,000,000 total reads/1,000,000 = 150 rPM
-OQ818329/RR034B_290_NODE_2      20,903,666 total reads/1,000,000 = 20.903 rPM

Partial sequences (accession number and CZID ID): 
-OQ818330/RR034B_002_NODE_4      81,508,414 total reads/1,000,000 = 81.508 rPM
-OQ818331/RR034B_014_NODE_22     81,944,496 total reads/1,000,000 = 81.944 rPM
-OQ818332/RR034B_016_NODE_39     71,534,446 total reads/1,000,000 = 71.534 rPM
-OQ818333/RR034B_018_NODE_5      55,495,770 total reads/1,000,000 = 55.495 rPM
-OQ818334/RR034B_028_NODE_1050   59,486,786 total reads/1,000,000 = 59.486 rPM
-OQ818335/RR034B_029_NODE_116    62,438,032 total reads/1,000,000 = 62.438 rPM
-OQ818336/RR034B_308_NODE_1      8,324,112 total reads/1,000,000 = 8.324 rPM
-OQ818337/RR034B_079_NODE_4      117,276,390 total reads/1,000,000 = 117.276 rPM
-OQ818338/RR034B_338_NODE_17     14,419,932 total reads/1,000,000 = 14.419 rPM
-OQ818339/RR034B_094_NODE_56     32,858,460 total reads/1,000,000 = 32.858 rPM
-OQ818340/RR034B_137_NODE_352    88,959,312 total reads/1,000,000 = 88.959 rPM
-OQ818341/RR034B_155_NODE_2      33,448,246 total reads/1,000,000 = 33.448 rPM
-OQ818342/RR034B_190_NODE_884    39,315,598 total reads/1,000,000 = 39.315 rPM
-OQ818343/RR034B_393_NODE_55     16,151,594 total reads/1,000,000 = 16.151 rPM
-OQ818344/RR034B_412_NODE_6      15,266,680 total reads/1,000,000 = 15.266 rPM
-OQ818345/RR034B_233_NODE_1161   27,847,184 total reads/1,000,000 = 27.847 rPM
-OQ818346/RR034B_248_NODE_215    16,855,124 total reads/1,000,000 = 16.855 rPM
-Q818347/RR034B_271_NODE_34      22,420,522 total reads/1,000,000 = 22.420 rPM
-OQ818348/RR034B_481_NODE_12     22,122,724 total reads/1,000,000 = 22.122 rPM

Use the R code [here]() to generate coverage plots
