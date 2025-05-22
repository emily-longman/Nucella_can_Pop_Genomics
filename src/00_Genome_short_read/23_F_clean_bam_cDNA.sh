#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=clean_bam_cDNA

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 # on one node
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss 
#SBATCH --time=30:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=60G

# Request CPU
#SBATCH --cpus-per-task=6

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# This script with sort and remove duplicates in the the bam file for cDNA, then index it.

#--------------------------------------------------------------------------------

# Load modules  
spack load samtools@1.10
PICARD=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard.jar
qualimap=/netfiles/nunezlab/Shared_Resources/Software/qualimap_v2.2.1/qualimap

#--------------------------------------------------------------------------------

# Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly

# This is the location where the reference genome and all its indexes are stored.
REFERENCE=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly/Base_Genome_Softmask/N.canaliculata_assembly.fasta.softmasked

#--------------------------------------------------------------------------------

# Define parameters
CPU=6
echo "using #CPUs ==" $CPU
JAVAMEM=30G # Java memory

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER_SCRATCH

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "cDNA_bam_qualimap" ]
then echo "Working cDNA_bam_qualimap folder exist"; echo "Let's move on."; date
else echo "Working cDNA_bam_qualimap folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/cDNA_bam_qualimap; date
fi

#--------------------------------------------------------------------------------

# Sort with picard
# Notice that once the file has been sorted it is added the "srt" suffix
java -Xmx$JAVAMEM -jar $PICARD SortSam \
I=$WORKING_FOLDER_SCRATCH/cDNA_bam/Nucella.cDNA.bam \
O=$WORKING_FOLDER_SCRATCH/cDNA_bam/Nucella.cDNA.srt.bam \
SO=coordinate \
VALIDATION_STRINGENCY=SILENT


# Remove duplicates with picard
# Notice that once the file has duplicates removed it is added the "rmdp" suffix
java -Xmx$JAVAMEM -jar $PICARD MarkDuplicates \
I=$WORKING_FOLDER_SCRATCH/cDNA_bam/Nucella.cDNA.srt.bam \
O=$WORKING_FOLDER_SCRATCH/cDNA_bam/Nucella.cDNA.srt.rmdp.bam \
M=$WORKING_FOLDER_SCRATCH/mapping_stats/Nucella.cDNA.dupstat.txt \
VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

# Index with samtools
samtools index $WORKING_FOLDER_SCRATCH/cDNA_bam/Nucella.cDNA.srt.rmdp.bam

#--------------------------------------------------------------------------------

# Run QC on the bam file

# Lets do QC on the bam file
$qualimap bamqc \
-bam $WORKING_FOLDER_SCRATCH/cDNA_bam/Nucella.cDNA.srt.rmdp.bam \
-outdir $WORKING_FOLDER_SCRATCH/cDNA_bam_qualimap/Qualimap_Nucella.cDNA \
--java-mem-size=$JAVAMEM