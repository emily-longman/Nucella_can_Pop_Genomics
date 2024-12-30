# Nucella_can_Pop_Genomics

## Project Goal

Characterize the population genomics of the channeled dogwhelk, Nucella canaliculata. 

Analyze the fine scale demography of Nucella canaliculata across the species range using neutral markers to estimate the level of gene flow and migration between populations. Identify loci under selection by comparing genotypes with environmental variables using a seascape genomics approach.

## File Structure

The files in this project are organized in the following structure:
 - data/
     - raw/
     - processed/
 - src/
    - 01_fastq_to_bam/
    - 02_enviro/
 - output/
     - figures/
     - tables/
 - docs/
 - scratch/


## Part 1 - Fastq to bam

Download the raw fastq reads. 

#### 01 - Reads_qc

These scripts will check the quality of the raw reads.

01_fastqc.sh - Run the program fastQC to produce quality reports on each individual raw read.  

02_multiqc.sh - Run the program multiQC on the fastQC outputs to produce one quality report for all raw reads.

#### 02 - Trim_and_Map

01_trim_read.sh - Use the program fastp to trim adapters as well as trim the reads based on quality.

02_index_reference.sh - Index the masked reference genome using the program bwa mem2. 

03_map_reads.sh - Map the reads to the masked reference genome using the program bwa mem2.

04_clean_bams.sh - Clean the bam files by filtering, sorting, and removing duplicates using the programs Picard and samtools. Then index the bams with samtools. Additionally, produce quality reports for each bam file using the program qualimap. 

05_merge_bams.sh - Merge the bam files for each population across the two lanes of sequencing using the program samtools. Additionally, produce quality reports for each bam file using the program qualimap.

06_coverage.sh - Generate a summary report using qualimap that shows the coverage for each bam file.
