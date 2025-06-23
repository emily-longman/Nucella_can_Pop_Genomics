# Geographic Divergence in Population Genomics and Shell Morphology Reveal History of Glacial Refugia in a Coastal Dogwhelk

## Project Summary

We assembled the draft genome for the channeled dogwhelk, *Nucella canaliculata*, using both Oxford Nanopore Technologies (ONT) long reads and AVITI PE300 short reads, and analyzed the fine scale demography of the species using pooled-sequencing of 19 populations distributed along ~1,500 km of the west coast of North America. Additionally, we performed geometric morphometrics of shell morphology to determine if spatial patterns of morphology are concordant with the phylogeographic patterns.

All sequencing was performed at [DNA Technologies and Expression Analysis Core](https://dnatech.ucdavis.edu/) at the UC Davis Genome Center and the bioinformatics pipeline was completed on the Vermont Advanced Computing Center ([VACC](https://www.uvm.edu/vacc)).

### Research Questions

1) How is genetic diversity in *N. canaliculata* spatially structured along ~1,500 km of the west coast of North America? 
2) Do patterns of genetic diversity appear to be structured more strongly by contemporary or historical processes? 
3) Does variation in shell morphology suggest a similar biogeographic pattern of divergence as the genomic data? 
4) Lastly, what are the molecular bases for intraspecific variation in shell morphology? 

## File Structure

The files in this project are organized in the following structure. All files for arrays and metadata are in "guide_files". All subsequent directories will be generated therein. Raw pool-seq data for the population genomic analyses can be found under the SRA BioProject accession number PRJNA1276871.
 - data/
     - raw/
     - processed/
 - guide_files
 - src/
    - 01_assemble_genome/
    - 02_fastq_to_vcf/
    - 03_population_structure/
    - 04_morphology/
 - output/
     - figures/
     - tables/

## Part 1 - Assemble the Draft Genome

The *N. canaliculata* draft genome was based on DNA extracted from one adult female *N. canaliculata* collected in June 2022 from Bodega Marine Reserve, California, USA. We used a hybrid assembly approach to assemble the genome, utilizing the continuinty of the Oxford Nanopore Technologies (ONT) PromethION long reads and the quality of the AVITI PE300 short reads. In order to optimze the assembly, in many of the steps, we generate multiple assemblies and then compare them. We determine the assembly quality and completeness we used [Quast](https://github.com/ablab/quast) and [BUSCO](https://www.expasy.org/resources/busco). 

### 01 - Prepare the Raw Data 

These scripts will format and filter the raw Oxford Nanopore Technologies (ONT) data and will trim the AVITI PE300 short reads. 

01_cat_reads.sh - Concatenate Oxford Nanopore Technologies (ONT) reads from 5 PromethION flow cells.   

02_filter_ONT_array.sh - Use the program [FiltLong](https://github.com/rrwick/Filtlong) to filter the ONT data by length. 

03_trim_reads.sh - Clean the AVITI short reads using [fastp](https://github.com/OpenGene/fastp).   

### 02 - DBG2OLC

These scripts will use a hybrid assembly approach using DBG2OLC to assemble a N. canaliculata genome. This approach takes advantage of the continuity of the ONT reads and the quality of the AVITI short reads.

01_sparseassembler_array.sh - Use the AVITI short reads to assemble contigs with [SparseAssembler](https://github.com/yechengxi/SparseAssembler).

02_DBG2OLC_array.sh - Use the best assembly from SparseAssembler as well as the 2,000 filtered ONT data as the inputs to [DBG20LC](https://github.com/yechengxi/DBG2OLC).

03_consensus: These steps will perform the consensus steps of DBG2OLC. In order to run effectively on the VACC, these steps were perform on partitions/chunks of the assembly. These steps also require additional consensus scripts that are in the folder "consensus_scripts" that is within '03_consensus'. These scripts are modified versions of split_and_run_sparc, split_reads_by_backbone.py and SeqIO.py which are provided by DBG2OLC.

- 01_convert_ONT_FQtoFA.sh - Convert the filtered ONT data from fastq to fasta format. 

- 02_cat_contigs_ONT.sh - Concatenate the SparseAssembler contigs and filtered ONT reads.

- 03_partition_guide_file.R - Create a partition file which will be used in the next step to break the contigs file and backbone into 50 contig chunks. 

- 04_partition_genome.sh - Partition the contigs file and backbone into 50 contig chunks.

- 05_guide_file_pt1 - (interactive session) Create a file with the backbone names.

- 05_guide_file_pt2.sh - Create array guide file for consensus step. 

- 06_run_consensus_caller.pt.sh - Run the first step in the consensus script (i.e., run the accompanying script "split_and_run_sparc.pt1")

- 06_run_consensus_caller.pt2_array.sh - Run the second step in the consensus script (i.e., run the accompanying script "split_and_run_sparc.pt2_array.sh")

- 06_run_consensus_caller.pt3.sh - This script will run the third step in the consensus script.

- 07_Quast.sh - This script will run Quast on the final consensus assembly.

- 08_BUSCO.sh - This script runs on BUSCO on the final assembly.

### 03 - Scaffold with ntlink

These scripts will use [ntLink](https://github.com/bcgsc/ntLink) to scaffold the hybrid assembly then assess the quality and completeness using Quast and BUSCO.

01_ntlink_array.sh - This script will use [ntLink](https://github.com/bcgsc/ntLink) to scaffold the assembly.

02_Quast_array.sh - This script will run Quast on the assembly generated from ntlink in the previous script.

03_BUSCO.sh - This script runs on BUSCO on the best ntlink assembly.

### 04 - Polish with Pilon

These scripts will use [Pilon](https://github.com/broadinstitute/pilon) to polish the genome 5 times with the AVITI reads. To be computationally efficient, the genome and bam files were broken into individual scaffolds and then polished individually using an array and loop structure. Each round, it will index the genome, map the short reads to the genome from the previous iteration, clean the bam files, generate a guide file with the scaffold names for the array, polish each scaffold, then concatenate the genome back together. The same 6 steps are peformed each round of polishing, thus they are only listed below once, rather than iterated 6 times.

01_index_genome.sh - Index the genome that was assembled by ntlink or the previous round of polishing.

02_map_reads.sh - Map the AVITI short reads to the genome.

03_clean_bams.sh - Clean the bam file.

04_create_guide_file.pt1 - (Interactive session) This code will extract the scaffold names from the genome.

04_create_guide_file.pt2 - Create a guide file of the scaffold names and group them into 30 scaffold chunks.

05_pilon_polish_chunks.sh - Polish the genome with Pilon using a combination of an array and a while loop.

06_cat_scaffolds.sh - Concatenated the polished scaffolds together to create a final assembly.

### 05 - Finalize the Genome

These sets of scripts will rename the scaffolds so that they are shorter and more meaningful and then will run Quast and BUSCO on the final assembly.

01_rename_scaffolds_guide_file_pt1 - (Interactive session) This code will extract the scaffold names from the genome.

01_rename_scaffolds_guide_file_pt2.R - Create a guide file of the scaffold names and group them into 30 scaffold chunks.

02_rename_scaffolds.sh - Use an array and loop to rename the individual scaffolds.

03_cat_scaffolds.sh - Concatenate the scaffolds together to produce a final assembly.

04_Quast.sh - Run Quast on the final assembly.

05_BUSCO.sh - Run BUSCO on the final assembly.

### 06 - Mask Repeats

Identify repeats and low complexity DNA sequences and mask them (i.e., replace with Ns).

01_repeat_softmask_C.gigas.sh - Use [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker) to softmask the genome - trained on the Pacific oyster, *Crassostrea gigas*. 

### 07 - Annotate the Genome

Annotate the softmasked genome with [Augustus](https://github.com/Gaius-Augustus/Augustus) and generate a [SnpEff](https://pcingola.github.io/SnpEff/) database. Note: prior to these steps, the genome was moved to netfiles (long term storage), thus path designation for the genome have changed.

01_Augustus.sh - Use [Augustus](https://github.com/Gaius-Augustus/Augustus) to annotate the softmasked genome (note model is trained on Drosophila melanogaster).

02_snpeff.sh - Use [SnpEff](https://pcingola.github.io/SnpEff/) to generate a SnpEff database that can be used to annotate subsequent VCF files.

## Part 2 - Fastq to VCF

This set of scripts will take the raw Pool-seq reads (PE150 of two lanes of NovaSeqX 25B) trim them, map them to the softmasked draft genome assembled in Part 1, and then call variants. Note: Prior to these steps the genome was moved to long term file storage (netfiles), thus path designation for the genome has changed.

### 01 - QC Reads

These scripts will check the quality of the raw pool-seq reads.

01_fastqc.sh - Run the program [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to produce quality reports on each individual raw read.  

02_multiqc.sh - Run the program [multiQC](https://seqera.io/multiqc/) on the fastQC outputs to produce one quality report for all raw reads.

### 02 - Trim and Map Reads

These scripts will trim the raw reads using [fastp](https://github.com/OpenGene/fastp), then map them to the genome using [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2).

01_trim_read.sh - Use the program [fastp](https://github.com/OpenGene/fastp) to trim adapters as well as trim the reads based on quality.

02_index_reference.sh - Index the masked reference genome using the program [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2). 

03_map_reads.sh - Map the reads to the masked reference genome using the program [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2).

04_clean_bams.sh - Clean the bam files by filtering, sorting, and removing duplicates using the programs [Picard](https://broadinstitute.github.io/picard/) and [Samtools](https://www.htslib.org/). Then index the bams with samtools. Additionally, produce quality reports for each bam file using the program [QualiMap](http://qualimap.conesalab.org/). 

05_merge_bams.sh - Merge the bam files for each population across the two lanes of sequencing using the program [Samtools](https://www.htslib.org/). Additionally, produce quality reports for each bam file using the program [QualiMap](http://qualimap.conesalab.org/).

06_bam_qc.sh - Run [QualiMap](http://qualimap.conesalab.org/) on the final bams. 

07_coverage.sh - Generate a summary report using [QualiMap](http://qualimap.conesalab.org/) that shows the coverage for each bam file.

### 03 - Call Variants

These scripts will identify variants for the 19 bam files using [freebayes](https://github.com/freebayes/freebayes) and generate a genome-wide panel of genetic variation. 

01_prep_bams.sh - This script will prep the bam files. More specifically, it will add read group information and index the final bam files. 

02_create_guide_file.pt1 - (Interactive session) This code will extract the scaffold names from the genome.

02_create_guide_file.pt2.R - Create a guide file of the scaffold names and group them into 30 scaffold chunks.

03_variant_calling.sh - This script will use [freebayes](https://github.com/freebayes/freebayes) to call variants. Note, chunk 617 required more than the allotted 30 hours and thus needed to be run again on a different partition on the VACC.

04_cat_vcf.sh - This script will concatenate together the chunked VCF files generated in the previous step from [freebayes](https://github.com/freebayes/freebayes).

05_pop_names.sh - This script will generate a guide file with the population names in the same order as the VCF.

06_filter_vcf.sh - This script will filter the VCF file created in the previous script. 

07_create_GDS.R - Create GDS object from the VCF.

## Part 3 - Population structure and genomic diversity

01_site_map.R - Generate site map for the 19 field sites distributed along ~1,500 km of the west coast of North America.

02_poolfstat.R - Use [poolfstat](https://cran.r-project.org/web/packages/poolfstat/index.html) to filter the SNP list and covert the VCF file to a pooldata object. Then, calculate FST statistics (a global FST metric, a block-jacknife estimation of FST, a multi-locus FST over a sliding window, and pair-wise population FST) and heterozygosity and generate a PCA of all SNPs. 

03_baypass - Use [BayPass](https://forge.inrae.fr/mathieu.gautier/baypass_public) to characterize population demography by analyzing the omega matrix.
- 01_format_baypass.R - Convert pooldata object to BayPass input format. 
- 02_generate_omega.sh - Generate an omega matrix using BayPass.
- 03_analyze_omega.R - Analyze the omega matrix and graph as a hierarchical clustering tree.

04_IBD.R - Test for isolation by distance (IBD) using Mantel tests.

05_npstat - Use [npstat](https://github.com/lucaferretti/npstat) to calculate nucleotide diversity, Watterson's estimator, and Tajima's D for each *N. canaliculata* population on a sliding window.  
- 01_create_guide_file.pt1 - (Interactive session) This code will extract the scaffold names from the genome.
- 01_create_guide_file.pt2.R -  Create a guide file of the scaffold names and group them into 50 scaffold chunks.
- 02_npstat.sh - Calculate diversity statistics (nucleotide diversity, Watterson's estimator, and Tajima's D) using npstat on a sliding window (length 25kb)
- 03_npstat_merge.sh - Merge the npstat diversity statistics txt files together.
- 04_npstat.R - Analyze and graph the diversity statistics.

## Part 4 - Analyses of shell morphology

Intraspecific variation in shell morphology was analyzed using landmark analysis based on photographs of the ventral surface of the dogwhelks. 15 landmarks were placed on the images using [tpsUtil](https://www.sbmorphometrics.org/soft-utility.html) and [tpsDig](https://www.sbmorphometrics.org/soft-dataacq.html). Some of the subsequent analyses were performed using the program [MorphoJ](https://morphometrics.uk/MorphoJ_page.html). All remaining geometric morphometric analyses and graphing were performed below in R. Additionally, we used [BayPass](https://forge.inrae.fr/mathieu.gautier/baypass_public) to identify outlier SNPs associated with variation in shell morphology.

01_geometric_morphometrics
- 01_morphology.R - Perform a Procrustes transformation on the landmarks then analyze and graph the morphology data with [geomorph](https://cran.r-project.org/web/packages/geomorph/index.html). Additionally, perform a Procrustes ANOVA with 1,000 permutationt o analyze the effects of size (log of the centroid size), site or latitude, and their interaction on dogwhelk shell shape.

02_baypass - Use [BayPass](https://forge.inrae.fr/mathieu.gautier/baypass_public) to identify outlier SNPs associated with morphological variation.
- 01_baypass_morphlogy_CV1.sh - Run baypass using the coefficient of variation of the PC1 of the morphology data as phenotypic input.
- 01_baypass_morphlogy_CV2.sh - Run baypass using the coefficient of variation of the PC2 of the morphology data as phenotypic input.
- 02_baypass_morphlogy_CV1_and_CV2.R - Graph baypass outputs.
- 03_snpeff_vcf.sh - Annotate VCF file using [SnpEff](https://pcingola.github.io/SnpEff/) database generated in Part 1.
- 04_annotate_SNPs.R - Annotate BayPass SNP lists.
- 05_filter.R - Filter the pooldata object for only the SNPs within the protein of interest.
- 06_top_SNPs.R - Analyses of the allele frequencies of the candidate loci.