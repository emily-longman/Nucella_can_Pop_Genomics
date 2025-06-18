#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=variant_calling_freebayes

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=30:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=60G 

# Request CPU
#SBATCH --cpus-per-task=4

# Submit job array
#SBATCH --array=1-631%25

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will use freebayes to call variants.

# Note: Since the N. canaliculata genome consistent of many scaffolds, I chunked the genome into 50 scaffolds partitions to utilize an array structure.
# For each slurm array ID, freebayes is running over the 50 scaffolds using a loop then the vcf files are combined together.  

#--------------------------------------------------------------------------------

# Load modules 
module load gcc/13.3.0-xp3epyt
module load freebayes/1.3.6-r67va2b
module load bcftools/1.19-iq5mwek

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

# This is the location where the reference genome is stored.
REFERENCE_FOLDER=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Oct2024/Crassostrea_softmask

# This is the path to the reference genome.
REFERENCE=$REFERENCE_FOLDER/N.canaliculata_assembly.fasta.softmasked.fa

#--------------------------------------------------------------------------------

## Import master partition file 
GUIDE_FILE=$WORKING_FOLDER/guide_files/Scaffold_names_array.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers. (dimensions: 2, 18919; 631 partitions)
# Scaffold name       # Partition
# Backbone_10001            1
# Backbone_10003            1
# Backbone_10004            1
# ....

#--------------------------------------------------------------------------------

# Determine partition to process 

# Echo slurm array task ID
echo ${SLURM_ARRAY_TASK_ID}

# Change directory 
cd $WORKING_FOLDER/data/processed/fastq_to_vcf

# Using the guide file, extract the scaffold names associated based on the Slurm array task ID for a given partition
awk '$2=='${SLURM_ARRAY_TASK_ID}'' $GUIDE_FILE | awk '{print $1}' > genome.scaffold.names.${SLURM_ARRAY_TASK_ID}.txt

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/fastq_to_vcf

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "vcf_freebayes" ]
then echo "Working vcf_freebayes folder exist"; echo "Let's move on."; date
else echo "Working vcf_freebayes folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/fastq_to_vcf/vcf_freebayes; date
fi

cd $WORKING_FOLDER/data/processed/fastq_to_vcf/vcf_freebayes

if [ -d "partitions" ]
then echo "Working partitions folder exist"; echo "Let's move on."; date
else echo "Working partitions folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/fastq_to_vcf/vcf_freebayes/partitions; date
fi

cd $WORKING_FOLDER/data/processed/fastq_to_vcf/vcf_freebayes/partitions

if [ -d "genome.scaffold.names.${SLURM_ARRAY_TASK_ID}" ]
then echo "Working genome.scaffold.names.${SLURM_ARRAY_TASK_ID} folder exist"; echo "Let's move on."; date
else echo "Working genome.scaffold.names.${SLURM_ARRAY_TASK_ID} folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/fastq_to_vcf/vcf_freebayes/partitions/genome.scaffold.names.${SLURM_ARRAY_TASK_ID}; date
fi

#--------------------------------------------------------------------------------

# Prepare bamlist
# This is a file with the name and full path of all the bam files to be processed.

# Move to bams folder
cd $WORKING_FOLDER/data/processed/fastq_to_vcf/RGSM_final_bams

# Create bamlist for all Nucella samples
ls *.bam > $WORKING_FOLDER/data/processed/fastq_to_vcf/guide_files/Nucella_pops_bam_names.list

BAMLIST=$WORKING_FOLDER/data/processed/fastq_to_vcf/guide_files/Nucella_pops_bam_names.list

#--------------------------------------------------------------------------------

# Use freebayes to call variants

# Note: 
# The freebayes algorithm exploits a neutral model of allele diffusion to impute most-confident genotypings across the
# entire population. In practice, the discriminant power of the method will improve if you run multiple samples simultaneously. 
# In other words, if your study has multiple individuals, you should run freebayes against them at the same time via a bamlist. 

# For the scaffolds in a given partition, generate a vcf file for each scaffold then cat them together. 

# Cat file of scaffold names and start while loop
cat $WORKING_FOLDER/data/processed/fastq_to_vcf/genome.scaffold.names.${SLURM_ARRAY_TASK_ID}.txt | \
while read SCAFFOLD 
do echo ${SCAFFOLD}

freebayes -f $REFERENCE -L $BAMLIST -r ${SCAFFOLD} -K -F 0.01 -C 1 -G 5 --limit-coverage 500 -n 4 -m 30 -q 20 | gzip -c > $WORKING_FOLDER/data/processed/fastq_to_vcf/vcf_freebayes/partitions/genome.scaffold.names.${SLURM_ARRAY_TASK_ID}/Nucella.freebayes.${SCAFFOLD}.vcf.gz

done 

# -K --pooled-continuous : Output all alleles which pass input filters, regardles of genotyping outcome or model.
## Least expensive option in time/memory because if we use pool-discrete it is necessary to specify a haploid pool size and it calculates the likelihoods for each possible configuration which can take a long time.

# -F --min-alternate-fraction N : Require at least this fraction of observations supporting an alternate allele within a single individual in order to evaluate the position.  default: 0.05
## In the case of Pool-seq it is better to reduce this value so as not to penalize rare alleles present in a pool.
## With the default value, any variant represented by <5% of reads in the pool would be ignored. 

# -C --min-alternate-count N : Require at least this count of observations supporting an alternate allele within a single individual in order to evaluate the position.  default: 2
## Reduce to 1 for the same reason as -F

# -G --min-alternate-total N : Require at least this count of observations supporting an alternate allele within the TOTAL population in order to use the allele in analysis.  default: 1
## Increased to 5. Given that we generally apply to many pools, avoid carrying around ultra rare alleles, which are probably sequencing errors. This also reduces computational costs.
## Caution however may have an impact for methods based on SFS.

# --limit-coverage N: Downsample per-sample coverage to this level if greater than this coverage. default: no limit

# -m --min-mapping-quality Q : Exclude alignments from analysis if they have a mapping quality less than Q.  default: 1
# -q --min-base-quality Q : Exclude alleles from analysis if their supporting base quality is less than Q.  default: 0
## Classic quality criteria (Note: may need to be readjusted if we recalibrate the Base Quality BQSR)?

# -n --use-best-n-alleles N : Evaluate only the best N SNP alleles, ranked by sum of supporting quality scores.  (Set to 0 to use all (i.e., default))
## Reduced to 4 (see recommendation for computational cost reduction). Be careful however because the other alleles are simply not considered (i.e., the real DP may vary from the actual DP, i.e., sum of coverage of the reported alleles)
## Most of the time, we will only retain bi-allelics (or very poorly covered alleles will be removed: see poolfstat)(in GATK HaplotypeCaller it is 7 by default because --max-alternate-alleles 6)

# --strict-vcf Generate strict VCF format (FORMAT/GQ will be an int)
# Genotype quality (not useful here since pooled-continuous) are given by default in log10 scale likelihood (=> heavy the output file)

#--------------------------------------------------------------------------------

# Combine the vcf files from all of the scaffolds in a given partition

bcftools concat \
$WORKING_FOLDER/data/processed/fastq_to_vcf/vcf_freebayes/partitions/genome.scaffold.names.${SLURM_ARRAY_TASK_ID}/*.vcf.gz \
-Oz -o $WORKING_FOLDER/data/processed/fastq_to_vcf/vcf_freebayes/partitions/genome.scaffold.names.${SLURM_ARRAY_TASK_ID}.vcf.gz

#--------------------------------------------------------------------------------

# Change directory
cd $WORKING_FOLDER/data/processed/fastq_to_vcf

# Housekeeping
rm genome.scaffold.names.${SLURM_ARRAY_TASK_ID}.txt
rm -R $WORKING_FOLDER/data/processed/fastq_to_vcf/vcf_freebayes/partitions/genome.scaffold.names.${SLURM_ARRAY_TASK_ID}

#--------------------------------------------------------------------------------

date
echo "done"

#--------------------------------------------------------------------------------

#### INFO and FORMAT field of the vcf

##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">
##INFO=<ID=DPB,Number=1,Type=Float,Description="Total read depth per bp at the locus; bases in reads overlapping / bases in haplotype">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">
##INFO=<ID=RO,Number=1,Type=Integer,Description="Count of full observations of the reference haplotype.">
##INFO=<ID=AO,Number=A,Type=Integer,Description="Count of full observations of this alternate haplotype.">
##INFO=<ID=PRO,Number=1,Type=Float,Description="Reference allele observation count, with partial observations recorded fractionally">
##INFO=<ID=PAO,Number=A,Type=Float,Description="Alternate allele observations, with partial observations recorded fractionally">
##INFO=<ID=QR,Number=1,Type=Integer,Description="Reference allele quality sum in phred">
##INFO=<ID=QA,Number=A,Type=Integer,Description="Alternate allele quality sum in phred">
##INFO=<ID=PQR,Number=1,Type=Float,Description="Reference allele quality sum in phred for partial observations">
##INFO=<ID=PQA,Number=A,Type=Float,Description="Alternate allele quality sum in phred for partial observations">
##INFO=<ID=SRF,Number=1,Type=Integer,Description="Number of reference observations on the forward strand">
##INFO=<ID=SRR,Number=1,Type=Integer,Description="Number of reference observations on the reverse strand">
##INFO=<ID=SAF,Number=A,Type=Integer,Description="Number of alternate observations on the forward strand">
##INFO=<ID=SAR,Number=A,Type=Integer,Description="Number of alternate observations on the reverse strand">
##INFO=<ID=SRP,Number=1,Type=Float,Description="Strand balance probability for the reference allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SRF and SRR given E(SRF/SRR) ~ 0.5, derived using Hoeffding's inequality">
##INFO=<ID=SAP,Number=A,Type=Float,Description="Strand balance probability for the alternate allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SAF and SAR given E(SAF/SAR) ~ 0.5, derived using Hoeffding's inequality">
##INFO=<ID=AB,Number=A,Type=Float,Description="Allele balance at heterozygous sites: a number between 0 and 1 representing the ratio of reads showing the reference allele to all reads, considering only reads from individuals called as heterozygous">
##INFO=<ID=ABP,Number=A,Type=Float,Description="Allele balance probability at heterozygous sites: Phred-scaled upper-bounds estimate of the probability of observing the deviation between ABR and ABA given E(ABR/ABA) ~ 0.5, derived using Hoeffding's inequality">
##INFO=<ID=RUN,Number=A,Type=Integer,Description="Run length: the number of consecutive repeats of the alternate allele in the reference genome">
##INFO=<ID=RPP,Number=A,Type=Float,Description="Read Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between RPL and RPR given E(RPL/RPR) ~ 0.5, derived using Hoeffding's inequality">
##INFO=<ID=RPPR,Number=1,Type=Float,Description="Read Placement Probability for reference observations: Phred-scaled upper-bounds estimate of the probability of observing the deviation between RPL and RPR given E(RPL/RPR) ~ 0.5, derived using Hoeffding's inequality">
##INFO=<ID=RPL,Number=A,Type=Float,Description="Reads Placed Left: number of reads supporting the alternate balanced to the left (5') of the alternate allele">
##INFO=<ID=RPR,Number=A,Type=Float,Description="Reads Placed Right: number of reads supporting the alternate balanced to the right (3') of the alternate allele">
##INFO=<ID=EPP,Number=A,Type=Float,Description="End Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between EL and ER given E(EL/ER) ~ 0.5, derived using Hoeffding's inequality">
##INFO=<ID=EPPR,Number=1,Type=Float,Description="End Placement Probability for reference observations: Phred-scaled upper-bounds estimate of the probability of observing the deviation between EL and ER given E(EL/ER) ~ 0.5, derived using Hoeffding's inequality">
##INFO=<ID=DPRA,Number=A,Type=Float,Description="Alternate allele depth ratio.  Ratio between depth in samples with each called alternate allele and those without.">
##INFO=<ID=ODDS,Number=1,Type=Float,Description="The log odds ratio of the best genotype combination to the second-best.">
##INFO=<ID=GTI,Number=1,Type=Integer,Description="Number of genotyping iterations required to reach convergence or bailout.">
##INFO=<ID=TYPE,Number=A,Type=String,Description="The type of allele, either snp, mnp, ins, del, or complex.">
##INFO=<ID=CIGAR,Number=A,Type=String,Description="The extended CIGAR representation of each alternate allele, with the exception that '=' is replaced by 'M' to ease VCF parsing.  Note that INDEL alleles do not have the first matched base (which is provided by default, per the spec) referred to by the CIGAR.">
##INFO=<ID=NUMALT,Number=1,Type=Integer,Description="Number of unique non-reference alleles in called genotypes at this position.">
##INFO=<ID=MEANALT,Number=A,Type=Float,Description="Mean number of unique non-reference allele observations per sample with the corresponding alternate alleles.">
##INFO=<ID=LEN,Number=A,Type=Integer,Description="allele length">
##INFO=<ID=MQM,Number=A,Type=Float,Description="Mean mapping quality of observed alternate alleles">
##INFO=<ID=MQMR,Number=1,Type=Float,Description="Mean mapping quality of observed reference alleles">
##INFO=<ID=PAIRED,Number=A,Type=Float,Description="Proportion of observed alternate alleles which are supported by properly paired read fragments">
##INFO=<ID=PAIREDR,Number=1,Type=Float,Description="Proportion of observed reference alleles which are supported by properly paired read fragments">
##INFO=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum depth in gVCF output block.">
##INFO=<ID=END,Number=1,Type=Integer,Description="Last position (inclusive) in gVCF output record.">
##INFO=<ID=technology.ILLUMINA,Number=A,Type=Float,Description="Fraction of observations supporting the alternate observed in reads from ILLUMINA">


##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Number of observation for each allele">
##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Reference allele observation count">
##FORMAT=<ID=QR,Number=1,Type=Integer,Description="Sum of quality of the reference observations">
##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count">
##FORMAT=<ID=QA,Number=A,Type=Integer,Description="Sum of quality of the alternate observations">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum depth in gVCF output block.">

