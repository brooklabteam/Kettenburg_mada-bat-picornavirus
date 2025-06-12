Mapping reads for full and partial sequences
---

This tutorial outlines methods for mapping the reads for contigs pulled from CZID in order to check coverage. I outline directions below.

---

At this point, I have a folder in Geneious with all of the hiqual sequences I pulled out of the command line. I made a duplicate of it and labeled it as read mapping so it's all contained in one folder.  

Next you should go to CZID [here]() and go to the respective folder for the sample you're looking for. Once you go into each folder, select a specific sample that needs a read mapped. 

After clicking on a sample, go to the upper right corner of the screen to the "Download button", click "that "View results folder". Scroll all the way to the bottom to "Step Nonhost Fastq: Filter original fastq/fasta input files to only contain non-host reads processed by CZ ID". Download both files. They should have .R1 and .R2 at the end of them. 

Once these are downloaded, drag the files to the read mapping folder in Geneious. Geneious will ask you to confirm the sequencing technology and will ask whether you want to pair the reads. Pair them and select Illumina as the sequencing technology. 

Now you have your non-host reads for a given sample. Select the non-host reads file and your sequence of interest and click align/assemble -> map to reference. Make sure your sequence is the sequence of interest for the reads to map to. I left all the default settings as is then mapped the reads. 


In the resulting file, it will show you the coverage per position in the genome. We're looking for decent coverage of at least two reads per site, ideally with overlapping positions to cover all areas of the genome. If the sequence does not have sufficient coverage, I make a file labeled "discarded" and put the sequence and the corresponding read report files in there. These sequences will not be further analyzed. On those that do have good coverage, click the graphs tab on the bottom right and select "export" to get a csv file of the coverage per position. From this, you can make one per sequence and use these to make graphs in R.

These values are "coverage" and we want the data in reads per million. We can get the total reads from CZID, per sample, go to the upper right corner and click "sample info". Then, click the tab that says pipelines and you'll see the total reads value. That reads value should be divided by the scaling factor of (total reads overall/1,000,000). The resulting value will be used in the csv sheet from Geneious. Divide every genome position value by this scaling factor to transform coverage to reads per million. 

Full genome sequences (accession number and CZID ID): 
-OQ818316/RR034B_029_NODE_84                62,438,032 total reads/1,000,000 = 62.438 rPM
-OQ818317/RR034B_079_NODE_3                 117,276,390 total reads/1,000,000 = 117.276 rPM
-OQ818318/RR034B_094_NODE_3                 32,858,460 total reads/1,000,000 = 32.858 rPM
-OQ818320/RR034B_163_NODE_204               46,343,504 total reads/1,000,000 = 46.343 rPM
-OQ818321/RR034B_189_NODE_1                 32,464,964 total reads/1,000,000 = 32.464 rPM
-OQ818322/RR034B_422_NODE_1                 17,718,526 total reads/1,000,000 = 17.718 rPM
-OQ818323/RR034B_244_NODE_40                95,323,590 total reads/1,000,000 = 95.323 rPM
-OQ818324/RR034B_259_NODE_68                18,588,814 total reads/1,000,000 = 18.588 rPM
-OQ818328/RR034B_288_NODE_100               150,000,000 total reads/1,000,000 = 150 rPM
-OQ818329/RR034B_290_NODE_2                 20,903,666 total reads/1,000,000 = 20.903 rPM
-PP766456/KEL86kidney_Plate_25_C8_NODE_1    46,323,502 total reads/1,000,000 = 46.323 rPM
-PP766459/KEL168_NODE_1                     53,993,066 total reads/1,000,000 = 53.993 rPM
-PP766469/MIZ405_Plate_32_A10_NODE_1        145,168,640 total reads/1,000,000 = 145.168 rPM
-/RR034B_242_fecRa_S67                        63,545,254 total reads/1,000,000 = 63.545 rPM

Partial sequences (accession number and CZID ID): 
-OQ818319/RR034B_096_NODE_2                   44,295,222 total reads/1,000,000 = 44.295 rPM
-OQ818325/RR034B_268_NODE_77                  31,773,930 total reads/1,000,000 = 31.773 rPM
-OQ818337/RR034B_079_NODE_4                   117,276,390 total reads/1,000,000 = 117.276 rPM
-OQ818340/RR034B_137_NODE_352                 88,959,312 total reads/1,000,000 = 88.959 rPM
-OQ818342/RR034B_190_NODE_884                 39,315,598 total reads/1,000,000 = 39.315 rPM
-OQ818343/RR034B_393_NODE_55                  16,151,594 total reads/1,000,000 = 16.151 rPM
-OQ818344/RR034B_412_NODE_6                   15,266,680 total reads/1,000,000 = 15.266 rPM
-OQ818345/RR034B_233_NODE_1161                27,847,184 total reads/1,000,000 = 27.847 rPM
-OQ818346/RR034B_248_NODE_215                 16,855,124 total reads/1,000,000 = 16.855 rPM
-OQ818347/RR034B_271_NODE_34                  22,420,522 total reads/1,000,000 = 22.420 rPM
-OQ818348/RR034B_481_NODE_12                  22,122,724 total reads/1,000,000 = 22.122 rPM
-PP766449/ANG22_Plate45_C6_NODE_21            72,382,862 total reads/1,000,000 = 72.382 rPM
-PP766450/ANG141_NODE_69                      36,103,456 total reads/1,000,000 = 36.103 rPM
-PP766451/ANGB110_Plate114_B9_NODE_12         43,982,106 total reads/1,000,000 = 43.982 rPM
-PP766452/ANGB110_Plate114_B9_NODE_58         43,982,106 total reads/1,000,000 = 43.982 rPM
-PP766453/ANGB110_Plate114_B9_NODE_320        43,982,106 total reads/1,000,000 = 43.982 rPM
-PP766454/ANGB121_Plate114_B7_NODE_306        43,982,106 total reads/1,000,000 = 43.982 rPM
-PP766455/KEL72_Plate45_D6_NODE_1             87,672,094 total reads/1,000,000 = 87.672 rpM
-PP766457/RR034B_079_NODE_4_and_80_assembled  117,276,390 total reads/1,000,000 = 117.276 rPM
-PP766458/RR034B_079_NODE_80                  117,276,390 total reads/1,000,000 = 117.276 rPM
-PP766460/RR034B_137_NODE_665                 88,959,312 total reads/1,000,000 = 88.959 rPM
-PP766461/RR034B_130_NODE_28                  33,895,162 total reads/1,000,000 = 33.895 rPM
-PP766462/KEL273_NODE_18                      27,464,336 total reads/1,000,000 = 27.464 rPM
-PP766463/KEL273_NODE_36                      27,464,336 total reads/1,000,000 = 27.464 rPM
-PP766464/RR034B_393_NODE_52                  16,151,594 total reads/1,000,000 = 16.151 rPM
-PP766465/KEL298_NODE_52                      15,266,680 total reads/1,000,000 = 15.266 rPM
-PP766466/RR034B_412_NODE_9                   15,266,680 total reads/1,000,000 = 15.266 rPM
-PP766467/KEL367_Plate_42_C2_NODE_1           13,746,818 total reads/1,000,000 = 13.746 rPM
-PP766468/RR034B_271_NODE_136                 22,420,522 total reads/1,000,000 = 22.420 rPM
-PP766470/MIZ405_Plate_32_A10_NODE_12         145,168,640 total reads/1,000,000 = 145.168 rPM
-PP766471/MIZ405_Plate_32_A10_NODE_20         145,168,640 total reads/1,000,000 = 145.168 rPM
-PP766472/MIZ405_Plate_32_A10_NODE_22         145,168,640 total reads/1,000,000 = 145.168 rPM
-PP766473/MIZ405_Plate_32_A10_NODE_44         145,168,640 total reads/1,000,000 = 145.168 rPM
-PP766474/MIZ405_Plate_32_A10_NODE_85         145,168,640 total reads/1,000,000 = 145.168 rPM
-PP766475/MIZ405_Plate_32_A10_NODE_105        145,168,640 total reads/1,000,000 = 145.168 rPM
-PP766476/MIZ405_Plate_32_A10_NODE_108        145,168,640 total reads/1,000,000 = 145.168 rPM
-PP766477/MIZ405_Plate_32_A10_NODE_129        145,168,640 total reads/1,000,000 = 145.168 rPM
-/RR034B_239_fecRa_S61                        53,902,938 total reads/1,000,000 = 53.902 rPM
-/RR034B_281_fecRa_S29                        43,437,908 total reads/1,000,000 = 43.437 rPM

Once you have exported all the csv coverage files and converted to rPM for plotting using the above information, you can concatentate all the files into one and use that in the R code to plot using the following code in terminal. 

```
cat *.csv > all_genome_coverage.csv

```

Use the R code [here]() to generate coverage plots
