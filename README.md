# Geographic Divergence in Population Genomics and Shell Morphology Reveal History of Glacial Refugia in a Coastal Dogwhelk

## Project Summary

Assembled the draft genome for the channeled dogwhelk, *Nucella canaliculata*, and analyzed the fine scale demography using pooled-sequencing of 19 populations distributed along ~1,500 km of the west coast of North America. Additionally, we performed geometric morphometrics of shell morphology to determine if spatial patterns of morphology are concordant with the phylogeographic patterns.

All sequencing was performed at DNA Technologies and Expression Analysis Core at the UC Davis Genome Center and the bioinformatics pipeline was completed on the Vermont Advanced Computing Center ([VACC](https://www.uvm.edu/vacc)).

### Research Questions

1) How is genetic diversity in *N. canaliculata* spatially structured along ~1,500km of the west coast of North America? 
2) Do patterns of genetic diversity appear to be structured more strongly by contemporary or historical processes? 
3) Does variation in shell morphology suggest a similar biogeographic pattern of divergence as the genomic data? 
4) Lastly, what are the molecular bases for intraspecific variation in shell morphology? 

## File Structure

The files in this project are organized in the following structure. All subsequent directories will be produced therein.
 - data/
     - raw/
     - processed/
 - src/
    - 01_assemble_genome/
    - 02_fastq_to_vcf/
    - 03_population_structure/
    - 04_morphology/
 - guide_files
 - output/
     - figures/
     - tables/
 - scratch/

## Part 1 - Assemble the Draft Genome

The *N. canaliculata* draft genome was based on DNA extracted from one adult female *N. canaliculata* collected in June 2022 from Bodega Marine Reserve, California, USA.

#### 01 - Prepare the Raw Data 

These scripts will format and filter the raw Oxford Nanopore Technologies (ONT) data and will trim the AVITI short reads. 

01_cat_reads.sh - Concatenate Oxford Nanopore Technologies (ONT) reads from 5 PromethION flow cells.   

02_filter_ONT_array.sh - Use the program Filtong (https://github.com/rrwick/Filtlong) to filter the ONT data by length. 

03_trim_reads.sh - Clean the AVITI short reads using fastp (https://github.com/OpenGene/fastp).   

#### 02 - DBG2OLC

These scripts will use a hybrid assembly approach using DBG2OLC taking advantage of the continuity of the ONT reads and the quality of the AVITI short reads.

01_sparseassembler_array.sh - This script will use the AVITI short reads to assemble contigs with SparseAssembler (https://github.com/yechengxi/SparseAssembler).

02_DBG2OLC_array.sh - This script will use the best assembly from SparseAssembler as well as the 2000 filt ONT data as the inputs to DBG20LC (https://github.com/yechengxi/DBG2OLC).

03_consensus: These steps will perform the consensus steps of DBG2OLC. In order to run effectively on the VACC, it will perform these steps on partitions/chunks of the assembly. 
- 01_convert_ONT_FQtoFA.sh - This script will convert the filtered ONT data from fastq to fasta format. 

- 02_cat_contigs_ONT.sh - This script will concatenate the SparseAssembler contigs and filtered ONT reads.

- 03_partition_guide_file.R - Create partition file which will be used in the next step to break the contigs file and backbone into 50 contig chunks. 

- 04_partition_genome.sh - This script will partition the contigs file and backbone into 50 contig chunks.

- 05_guide_file_pt1 - (interactive session) This script will create a file with the backbone names.

- 05_guide_file_pt2.sh - Create array guide file for consensus step. 

- 06_run_consensus_caller.pt.sh - This script will run the first step in the consensus script. Specifically, it will run the consensus script "split_and_run_sparc.pt1".

- 06_run_consensus_caller.pt2_array.sh - This script will run the second step in the consensus script. Specifically it will run the consensus script "split_and_run_sparc.pt2_array.sh".

- 06_run_consensus_caller.pt3.sh - This script will run the third step in the consensus script.

- 07_Quast.sh - This script will run Quast on the final consensus assembly.

- 08_BUSCO.sh - This script runs on BUSCO on an assembly.

#### 03 - Scaffold with ntlink

These scripts will use a ntLink (https://github.com/bcgsc/ntLink) to scaffold the hybrid assembly then assess the quality and completeness using Quast and BUSCO.

01_ntlink_array.sh - This script will use ntLink (https://github.com/bcgsc/ntLink) to scaffold the assembly

02_Quast_array.sh - This script will run Quast on the assembly generated from ntlink in the previous script.

03_BUSCO.sh - This script runs on BUSCO on the best ntlink assembly.

#### 04 - Polish with Pilon

These scripts will use Pilon [v1.24] (https://github.com/broadinstitute/pilon) to polish the genome 5 times with the AVITI reads. To be computationally efficient, the genome and bam files were broken into individual scaffolds and then polished individually using an array and loop structure. Each round, it will index the genome, map the short reads to the genome from the previous iteration, clean the bam files, generate a guide file with the scaffold names for the array, polish each scaffold, then concatenate the genome back together. The same 6 steps are peformed each round of polishing, thus they are only listed below once, rather than iterated 6 times.

01_index_genome.sh - Index the genome that was assembled by ntlink or the previous round of polishing.

02_map_reads.sh - Map the AVITI short reads to the genome.

03_clean_bams.sh - Clean the bam file.

04_create_guide_file.pt1 - (Interactive session) This code will extract the scaffold names from the genome.

04_create_guide_file.pt2 - Create a guide file of the scaffold names and group them into 30 scaffold chunks.

05_pilon_polish_chunks.sh - Polish the genome with Pilon using a combination of an array and while loop.

06_cat_scaffolds.sh - The polished scaffolds are concatenated together to create a final assembly.

#### 05 - Finalize the Genome

These sets of scripts will rename the scaffolds so that they are shorter and more meaningful and then will run Quast and BUSCO on the final assembly.

01_rename_scaffolds_guide_file_pt1 - (Interactive session) This code will extract the scaffold names from the genome.

01_rename_scaffolds_guide_file_pt2.R - Create a guide file of the scaffold names and group them into 30 scaffold chunks.

02_rename_scaffolds.sh - Use an array and loop to rename the individual scaffolds

03_cat_scaffolds.sh - Concatenate the scaffolds together to produce a final assembly.

04_Quast.sh - Run Quast on the final assembly.

05_BUSCO.sh - Run BUSCO on the final assembly.

#### 06 - Mask Repeats

01_repeat_softmask_C.gigas.sh - Use repeat masker (https://github.com/Dfam-consortium/RepeatMasker) to softmask the genome. 

#### 07 - Annotate the Genome


## Part 2 - Fastq to VCF

This set of scripts will take the Pool-seq raw reads, trim them then map them to the genome and then call variants. 

#### 01 - QC reads

These scripts will check the quality of the pool-seq raw reads.

01_fastqc.sh - Run the program fastQC to produce quality reports on each individual raw read.  

02_multiqc.sh - Run the program multiQC on the fastQC outputs to produce one quality report for all raw reads.

#### 02 - Trim and Map Reads

01_trim_read.sh - Use the program fastp to trim adapters as well as trim the reads based on quality.

02_index_reference.sh - Index the masked reference genome using the program bwa mem2. 

03_map_reads.sh - Map the reads to the masked reference genome using the program bwa mem2.

04_clean_bams.sh - Clean the bam files by filtering, sorting, and removing duplicates using the programs Picard and samtools. Then index the bams with samtools. Additionally, produce quality reports for each bam file using the program qualimap. 

05_merge_bams.sh - Merge the bam files for each population across the two lanes of sequencing using the program samtools. Additionally, produce quality reports for each bam file using the program qualimap.

06_bam_qc.sh - Run qualimap on the final bams. 

07_coverage.sh - Generate a summary report using qualimap that shows the coverage for each bam file.

08_prep_bams.sh - This script will prep the bam files. More specifically, it will add read group information and index the final bam files. 

09_make_chunks - Create guide file for array.

10_variant_calling.sh - This script will use freebayes (https://github.com/freebayes/freebayes) to call variants.

11_cat_vcf.sh - This script will cat together the chunked vcf files generated in the previous step from freebayes.

12_pop_names.sh - This script will generate a guide file with the population names in the same order as the vcf.

13_filter_vcf.sh - This script will filter the vcf file created in the previous script. 

14_create_GDS.R - Create GDS object from the vcf.

## Part 3 - Population structure and genomic diversity

01_site_map.R - Generate site map. 

02_poolfstat.R - Use poolfstat (https://cran.r-project.org/web/packages/poolfstat/index.html) to calculate Fst statistics and generate PCA of all SNPs. 

03_baypass - Use BayPass (https://forge.inrae.fr/mathieu.gautier/baypass_public) to characterize population dmoegraphy by analyzing the omega matrix.
- 01_format_baypass.R
- 02_generate_omega.sh
- 03_analyze_omega.R

04_IBD.R - Test for isolation by distance (IBD). 

05_npstat - Use npstat (https://github.com/lucaferretti/npstat) to calculate nucleotide diversity, Watterson's estimator, and Tajima's D for each *N. canaliculata* population on a sliding window.  
- 01_create_guide_file
- 02_npstat.sh
- 03_npstat_merge.sh
- 04_npstat.R


## Part 4 - Shell morphology

01_geometric_morphometrics
- 01_morphology.R

02_baypass