
Finding positive picorna samples after mNGS+CZID
=======

1. First, we downloaded all non-host contigs derived from fecal, urine, tissue, or saliva and mapping to ANY taxon from CZID in the following projects: 
 - RR034B1_feces (sequenced at CZID)
 - RR034B2_urine_wholeblood_novaseq (sequenced at CZID)
 - RR034B_throat_swab_raw_RNASeq_NovaSeq (sequenced at CZID)
 - Cara_Brook_225L55LT3 (sequenced at NIH for a different project)
 - Cara_Brook_22C2NFLT3 (sequenced at NIH for a different project)
 - CaraBrook_Data_Nov2023_tissue (sequenced at NIH for a different project)
 - Cara_Brook_HK2N3DRX2 (sequenced at NIH for a different project)
 - Cara_Brook_HK2KCDRX2 (sequenced at NIH for a different project)
This can be done in bulk, manually, on CZID.net in the top right-hand corner.

(Note that we did not take contigs from HeLa controls or water, and in the case of the 'RR034B2_urine_wholeblood_novaseq' project, we only looked at urine samples)

Note that when you download all the non-host contigs, it will produce a folder with a separate fasta file for each sample, which lists the contigs by node number but does not include the sample ID. Before joining all the contigs (nodes) together, you need to distinguish them by sample ID. I wrote Rscripts that parses this for each filetype (rename-fastas-feces, etc.). To rename your files and the headers within them, copy the appropriate Rscript for the tissue type into your appropriate downloads folder, cd into that folder on the command line, and simply type:

```

Rscript rename-fastas-feces.R
Rscript rename-fastas-urine.R
Rscript rename-fastas-tissue.R

```

Once the file finishes running, concatenate all the abbreviated filenames into a compiled file for downstream analyses. Note that these files are too big to include in the GitHub repo. All NIH files are from a different project.

```

cat *.fasta>rr034b1_feces_all_non_host_contigs.fasta
cat *.fasta>rr034b2_urine_all_non_host_contigs.fasta
cat *.fasta>nih_urine_jun2023_all_non_host_contigs.fasta
cat *.fasta>nih_urine_dec2023_all_non_host_contigs1.fasta
cat *.fasta>nih_urine_dec2023_all_non_host_contigs2.fasta
cat *.fasta>nih_tissue_dec2023_all_non_host_contigs.fasta

```
---

2. Now, because many of the contigs will overlap in sequences and slow our mapping down, we can deduplicate the compiled .fasta for each of the contigs for all tissues using the program [CD-HIT](http://weizhong-lab.ucsd.edu/cd-hit/). You'll need to [install this](https://github.com/weizhongli/cdhit/wiki/2.-Installation) on your home computer first. If using MacOS, I recommend using bioconda to do it--see [here](https://anaconda.org/bioconda/cd-hit). Once installed, you can deduplicate each of the contig files via the following command line script in the same folder as the downloaded contigs. This command :

```
cd-hit-est -i rr034b1_feces_all_non_host_contigs.fasta -c 0.95 -o rr034b1_feces_all_non_host_contigs_DEDUP.fasta -M 0
cd-hit-est -i rr034b2_urine_all_non_host_contigs.fasta -c 0.95 -o rr034b2_urine_all_non_host_contigs_DEDUP.fasta -M 0
cd-hit-est -i nih_urine_jun2023_all_non_host_contigs.fasta -c 0.95 -o nih_urine_jun2023_all_non_host_contigs_DEDUP.fasta -M 0
cd-hit-est -i nih_throat_jun2023_all_non_host_contigs.fasta -c 0.95 -o nih_throat_jun2023_all_non_host_contigs_DEDUP.fasta -M 0
cd-hit-est -i nih_urine_dec2023_all_non_host_contigs1.fasta -c 0.95 -o nih_urine_dec2023_all_non_host_contigs_DEDUP1.fasta -M 0
cd-hit-est -i nih_urine_dec2023_all_non_host_contigs2.fasta -c 0.95 -o nih_urine_dec2023_all_non_host_contigs_DEDUP2.fasta -M 0
cat *.fasta>nih_urine_dec2023_all_non_host_contigs_DEDUP.fasta
cd-hit-est -i nih_tissue_dec2023_all_non_host_contigs.fasta -c 0.95 -o nih_tissue_dec2023_all_non_host_contigs_DEDUP.fasta -M 0
cd-hit-est -i gce_ns_dec2023_all_non_host_contigs.fasta -c 0.95 -o gce_ns_dec2023_all_non_host_contigs_DEDUP.fasta -M 0

```
I ran these on my local computer, but they will run much quicker on midway. Again, these files are quite large and were used for multiple viral discovery projects so are not on GitHub.

---

3. Next,  after the contigs are deduped (or simulataneously as you are doing this), you can download all the (a) nucleotide and (b) protein full genome **reference sequences** under 

(i) Picornavirales taxid 464095

(ii) Unclassified picornavirales taxid 675074

from [NCBI virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/). 

Concatenate the files where there are multiple sub-groups. 

---

4. Now, you need to make sure that the command line version of NCBI-Blast is installed on your computer. See [here](https://www.ncbi.nlm.nih.gov/books/NBK569861/) for directions. If you are using Mac OS X, running the .dmg installer will probably give you the most success. To test if your installation worked, from anywhere in the command line, try typing 'blastn' and hitting enter. If you receive the following message, then you should be good to go:

```
BLAST query/options error: Either a BLAST database or subject sequence(s) must be specified
```

At this point I was still running commands on my local computer, it is faster to do on Midway so if you go this route you will need to install Blast to the computing cluster.

---

4. Once Blast is installed, you need to make a reference database from your downloaded picornavirus sequences. To do this, on the command line, you can "cd" into the folder where these are contained and use the following command on each of the fasta files produced in step 2 to build these into two different blast databases:

```
makeblastdb -in NCBI_picornavirales_nt.fasta -dbtype nucl -parse_seqids -out picornavirales_nt
makeblastdb -in NCBI_picornavirales_protein.fasta -dbtype prot -parse_seqids -out picornavirales_protein

```

Note that you may have to delete and retype the dashes above in your own command line run. This may not copy/paste easily. The commands above should generate a suite of files in the same folder that have the prefix specified after the "out" command.

---

5. Now that the reference sequence is in hand and the non-host contigs deduplicated, you can kick off a command line blast for each contig subset on the two above databases. I ran these using SLURM scripts on the Midway computing cluster using the example script:

```
#!/bin/bash
#SBATCH --job-name=blastn-feces-picornaviridae
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --time=36:00:00

blastn -word_size 10 -evalue 0.001 -query rr034b1_feces_all_non_host_contigs_DEDUP.fasta -db picornaviridae_nt -outfmt '6 qseqid nident pident length evalue bitscore sgi sacc stitle'  -max_target_seqs 10 -out picornaviridae_blast_feces_nt.txt

```
Basically, I am re-doing a more focused version of CZID to see if any other "hits" to viruses specifically shake out.

You will run four BLASTs in total: 2 nucleotide and 2 protein BLASTs, one for each of the deduplicated contigs above. First, run a "blastn"" alignment of deduplicated set of contigs with the virus_nt database, then run a "blastx" alignment of the deduplicated set of contigs with the virus_aa database. Scripts for both are listed below (make sure that you upload ALL the virus_aa and virus_nt files into the same folder for this to be able to run):

Picornavirales
```
blastn -word_size 10 -evalue 0.001 -query rr034b1_feces_all_non_host_contigs_DEDUP.fasta -db picornavirales_nt -outfmt '6 qseqid nident pident length evalue bitscore sgi sacc stitle'  -max_target_seqs 10 -out picornavirales_blast_feces_nt.txt

blastx -word_size 3 -evalue 0.001 -query rr034b1_feces_all_non_host_contigs_DEDUP.fasta -db picornavirales_protein -outfmt '6 qseqid nident pident length evalue bitscore sgi sacc stitle' -max_target_seqs 10 -out picornavirales_blast_feces_protein.txt

blastn -word_size 10 -evalue 0.001 -query rr034b2_urine_all_non_host_contigs_DEDUP.fasta -db picornavirales_nt -outfmt '6 qseqid nident pident length evalue bitscore sgi sacc stitle'  -max_target_seqs 10 -out picornavirales_blast_urine_nt.txt

blastx -word_size 3 -evalue 0.001 -query rr034b2_urine_all_non_host_contigs_DEDUP.fasta -db picornavirales_protein -outfmt '6 qseqid nident pident length evalue bitscore sgi sacc stitle' -max_target_seqs 10 -out picornavirales_blast_urine_protein.txt

blastn -word_size 10 -evalue 0.001 -query rr034b_throat_swab_all_non_host_contigs_DEDUP.fasta -db picornavirales_nt -outfmt '6 qseqid nident pident length evalue bitscore sgi sacc stitle'  -max_target_seqs 10 -out picornavirales_blast_throat_nt.txt

blastx -word_size 3 -evalue 0.001 -query rr034b_throat_swab_all_non_host_contigs_DEDUP.fasta -db picornavirales_protein -outfmt '6 qseqid nident pident length evalue bitscore sgi sacc stitle' -max_target_seqs 10 -out picornavirales_blast_throat_protein.txt

blastn -word_size 10 -evalue 0.001 -query nih_urine_jun2023_all_non_host_contigs_DEDUP.fasta -db picorna_nt -outfmt '6 qseqid nident pident length evalue bitscore sgi sacc stitle'  -max_target_seqs 10 -out picornavirales_blast_urine_nt.txt

blastx -word_size 3 -evalue 0.001 -query nih_urine_jun2023_all_non_host_contigs_DEDUP.fasta -db picorna_aa -outfmt '6 qseqid nident pident length evalue bitscore sgi sacc stitle' -max_target_seqs 10 -out picornavirales_blast_urine_protein.txt

blastn -word_size 10 -evalue 0.001 -query nih_urine_dec2023_all_non_host_contigs_DEDUP.fasta -db picorna_nt -outfmt '6 qseqid nident pident length evalue bitscore sgi sacc stitle'  -max_target_seqs 10 -out picornavirales_blast_urine_nt.txt

blastx -word_size 3 -evalue 0.001 -query nih_urine_dec2023_all_non_host_contigs_DEDUP.fasta -db picorna_aa -outfmt '6 qseqid nident pident length evalue bitscore sgi sacc stitle' -max_target_seqs 10 -out picornavirales_blast_urine_protein.txt

blastn -word_size 10 -evalue 0.001 -query nih_tissue_dec2023_all_non_host_contigs_DEDUP.fasta -db picorna_nt -outfmt '6 qseqid nident pident length evalue bitscore sgi sacc stitle'  -max_target_seqs 10 -out picornavirales_blast_tissue_nt.txt

blastx -word_size 3 -evalue 0.001 -query nih_tissue_dec2023_all_non_host_contigs_DEDUP.fasta -db picorna_aa -outfmt '6 qseqid nident pident length evalue bitscore sgi sacc stitle' -max_target_seqs 10 -out picornavirales_blast_tissue_protein.txt

```

Again, these can be run locally, but I used Midway to speed things up. The searches did not take so long (minutes!) after the deduplication step above for the urine. The fecal search was longer.

---

6.  After the blast finishes, you'll want to curate a bit to the high quality hits. After Amy's lead, I went ahead and parsed for alignments that show alignment length > 100 aa and bit score > 100. 

Here are the script for the nt and aa parses::

Picornavirales
```
cat picornavirales_blast_feces_nt.txt | awk -F\t '($4>99 && $5<0.00001) {print $1,$3, $4, $5, $8,$9}' > picornavirales_blast_feces_nt_results_100len5eval.txt

cat picornavirales_blast_feces_protein.txt | awk -F\t '($4>99 && $6>99) {print $1,$3, $4,$5,$8,$9}' > picornavirales_blast_feces_aa_results_100len100bit.txt


cat picornavirales_blast_urine_nt.txt | awk -F\t '($4>99 && $5<0.00001) {print $1,$3, $4, $5, $8,$9}' > picornavirales_blast_urine_nt_results_100len5eval.txt

cat picornavirales_blast_urine_protein.txt | awk -F\t '($4>99 && $6>99) {print $1,$3, $4,$5,$8,$9}' > picornavirales_blast_urine_aa_results_100len100bit.txt

cat picornavirales_blast_throat_nt.txt | awk -F\t '($4>99 && $5<0.00001) {print $1,$3, $4, $5, $8,$9}' > picornavirales_blast_throat_nt_results_100len5eval.txt

cat picornavirales_blast_throat_protein.txt | awk -F\t '($4>99 && $6>99) {print $1,$3, $4,$5,$8,$9}' > picornavirales_blast_throat_aa_results_100len100bit.txt

cat picornavirales_blast_urine_nt.txt | awk -F\t '($4>99 && $5<0.00001) {print $1,$3, $4, $5, $8,$9}' > picornavirales_blast_urine_nt_results_100len5eval.txt

cat picornavirales_blast_urine_protein.txt | awk -F\t '($4>99 && $6>99) {print $1,$3, $4,$5,$8,$9}' > picornavirales_blast_urine_aa_results_100len100bit.txt

cat picornavirales_blast_throat_nt.txt | awk -F\t '($4>99 && $5<0.00001) {print $1,$3, $4, $5, $8,$9}' > picornavirales_blast_throat_nt_results_100len5eval.txt

cat picornavirales_blast_throat_protein.txt | awk -F\t '($4>99 && $6>99) {print $1,$3, $4,$5,$8,$9}' > picornavirales_blast_throat_aa_results_100len100bit.txt


cat picornavirales_blast_tissue_nt.txt | awk -F\t '($4>99 && $5<0.00001) {print $1,$3, $4, $5, $8,$9}' > picornavirales_blast_tissue_nt_results_100len5eval.txt

cat picornavirales_blast_tissue_protein.txt | awk -F\t '($4>99 && $6>99) {print $1,$3, $4,$5,$8,$9}' > picornavirales_blast_tissue_aa_results_100len100bit.txt

```


---

7. Once the blast results have been sub-selected a bit, you can summarize them to link back the hits to the samples of interest. Within the same folder as your output, try the following script to save the unique contig IDs which align to virus using the following scripts and subsequent command to save the unique sample IDs for sample types:

Hiqual summaries 

Picornavirales
```
cat picornavirales_blast_feces_nt_results_100len5eval.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_contigs_feces_nt_hiqual.txt

cat picornavirales_unique_contigs_feces_nt_hiqual.txt | awk -F\_ '{print $1"_"$2}' | sort | uniq > picornavirales_unique_sampleID_feces_nt_hiqual.txt

cat picornavirales_blast_feces_aa_results_100len100bit.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_contigs_feces_aa_hiqual.txt

cat picornavirales_unique_contigs_feces_aa_hiqual.txt | awk -F\_ '{print $1"_"$2}' | sort | uniq > picornavirales_unique_sampleID_feces_aa_hiqual.txt



cat picornavirales_blast_urine_nt_results_100len5eval.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_contigs_urine_nt_hiqual.txt

cat picornavirales_unique_contigs_urine_nt_hiqual.txt | awk -F\_ '{print $1"_"$2}' | sort | uniq > picornavirales_unique_sampleID_urine_nt_hiqual.txt

cat picornavirales_blast_urine_aa_results_100len100bit.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_contigs_urine_aa_hiqual.txt

cat picornavirales_unique_contigs_urine_aa_hiqual.txt | awk -F\_ '{print $1"_"$2}' | sort | uniq > picornavirales_unique_sampleID_urine_aa_hiqual.txt



cat picornavirales_blast_throat_nt_results_100len5eval.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_contigs_throat_nt_hiqual.txt

cat picornavirales_unique_contigs_throat_nt_hiqual.txt | awk -F\_ '{print $1"_"$2}' | sort | uniq > picornavirales_unique_sampleID_throat_nt_hiqual.txt

cat picornavirales_blast_throat_aa_results_100len100bit.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_contigs_throat_aa_hiqual.txt

cat picornavirales_unique_contigs_throat_aa_hiqual.txt | awk -F\_ '{print $1"_"$2}' | sort | uniq > picornavirales_unique_sampleID_throat_aa_hiqual.txt



cat picornavirales_blast_urine_nt_results_100len5eval.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_contigs_urine_nt_hiqual.txt

cat picornavirales_unique_contigs_urine_nt_hiqual.txt | awk -F\_ '{print $1"_"$2}' | sort | uniq > picornavirales_unique_sampleID_urine_nt_hiqual.txt

cat picornavirales_blast_urine_aa_results_100len100bit.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_contigs_urine_aa_hiqual.txt

cat picornavirales_unique_contigs_urine_aa_hiqual.txt | awk -F\_ '{print $1"_"$2}' | sort | uniq > picornavirales_unique_sampleID_urine_aa_hiqual.txt



cat picornavirales_blast_throat_nt_results_100len5eval.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_contigs_throat_nt_hiqual.txt

cat picornavirales_unique_contigs_throat_nt_hiqual.txt | awk -F\_ '{print $1"_"$2}' | sort | uniq > picornavirales_unique_sampleID_throat_nt_hiqual.txt

cat picornavirales_blast_throat_aa_results_100len100bit.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_contigs_throat_aa_hiqual.txt

cat picornavirales_unique_contigs_throat_aa_hiqual.txt | awk -F\_ '{print $1"_"$2}' | sort | uniq > picornavirales_unique_sampleID_throat_aa_hiqual.txt



cat picornavirales_blast_tissue_nt_results_100len5eval.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_contigs_tissue_nt_hiqual.txt

cat picornavirales_unique_contigs_tissue_nt_hiqual.txt | awk -F\_ '{print $1"_"$2}' | sort | uniq > picornavirales_unique_sampleID_tissue_nt_hiqual.txt

cat picornavirales_blast_tissue_aa_results_100len100bit.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_contigs_tissue_aa_hiqual.txt

cat picornavirales_unique_contigs_tissue_aa_hiqual.txt | awk -F\_ '{print $1"_"$2}' | sort | uniq > picornavirales_unique_sampleID_tissue_aa_hiqual.txt

```

Broad catch-all summaries

Picornavirales
```
cat picornavirales_blast_feces_nt.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_contigs_feces_nt.txt

cat picornavirales_blast_feces_protein.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_contigs_feces_aa.txt

cat picornavirales_blast_urine_nt.txt | awk -F\_ '{print $1"_"$2}' | sort | uniq > picornavirales_unique_contigs_urine_nt.txt

cat picornavirales_blast_urine_protein.txt | awk -F\_ '{print $1"_"$2}' | sort | uniq > picornavirales_unique_contigs_urine_protein.txt

cat picornavirales_blast_throat_nt.txt | awk -F\_ '{print $1"_"$2}' | sort | uniq > picornavirales_unique_contigs_throat_nt.txt

cat picornavirales_blast_throat_protein.txt | awk -F\_ '{print $1"_"$2}' | sort | uniq > picornavirales_unique_contigs_throat_protein.txt



cat picornavirales_unique_contigs_feces_nt.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_sampleID_feces_nt.txt

cat picornavirales_unique_contigs_feces_aa.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_sampleID_feces_aa.txt

cat picornavirales_unique_contigs_urine_nt.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_sampleID_urine_nt.txt

cat picornavirales_unique_contigs_urine_protein.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_sampleID_urine_protein.txt

cat picornavirales_unique_contigs_throat_nt.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_sampleID_throat_nt.txt

cat picornavirales_unique_contigs_throat_protein.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_sampleID_throat_protein.txt


cat picornavirales_blast_urine_nt.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_contigs_urine_nt.txt

cat picornavirales_blast_urine_protein.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_contigs_urine_protein.txt

cat picornavirales_unique_contigs_urine_nt.txt | awk -F\_ '{print $1}' | sort | uniq > picornavirales_unique_sampleID_urine_nt.txt

cat picornavirales_unique_contigs_urine_protein.txt | awk -F\_ '{print $1}' | sort | uniq > picornavirales_unique_sampleID_urine_protein.txt



cat picornavirales_blast_throat_nt.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_contigs_throat_nt.txt

cat picornavirales_blast_throat_protein.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_contigs_throat_protein.txt

cat picornavirales_unique_contigs_throat_nt.txt | awk -F\_ '{print $1}' | sort | uniq > picornavirales_unique_sampleID_throat_nt.txt

cat picornavirales_unique_contigs_throat_protein.txt | awk -F\_ '{print $1}' | sort | uniq > picornavirales_unique_sampleID_throat_protein.txt



cat picornavirales_blast_tissue_nt.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_contigs_tissue_nt.txt

cat picornavirales_blast_tissue_protein.txt | awk '{print $1}' | sort | uniq > picornavirales_unique_contigs_tissue_protein.txt

cat picornavirales_unique_contigs_tissue_nt.txt | awk -F\_ '{print $1}' | sort | uniq > picornavirales_unique_sampleID_tissue_nt.txt

cat picornavirales_unique_contigs_tissue_protein.txt | awk -F\_ '{print $1}' | sort | uniq > picornavirales_unique_sampleID_tissue_protein.txt

```

---

8. For calling positives in cases where there was a discrepancy between this (stringent) pipeline and CZID, we will accept them as positive hits if (and only if!) the reads from that sample assembled into one or more contigs. In this case, contigs should only be acceptable if the average read depth at that contig is 2 or more reads (per Amy's rule). So, in manually curating any positive samples from CZID, check the broad (not hiqual) contig summary file for that sample and only call it as positive if it has at least one contig with >2 reads for average coverage.

Double checked every contig by using BLAST online for both nt and aa identities, then checked read depth through CZID by clicking on each sample ID, then downloading the R1 and R2 files at the very bottom in the results folder (they are called step nonhost Fastq), pairing, then mapping reads to the contig. 

Those contigs with low read support were not used for further analysis. I then separated the contigs into picornaviridae/caliciviridae and other picornavirales, then subset into full and partial sequences. All sample info and contigs sequences are found in the folder [here]().

