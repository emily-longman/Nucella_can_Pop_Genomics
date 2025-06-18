#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Move_cDNA

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss --7 day limit
#SBATCH --time=8:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=20G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# This code will cat all of the raw reads for the cDNA together 

Nucella_netfiles=/netfiles/pespenilab_share/Nucella/raw/

#--------------------------------------------------------------------------------

# Move to cDNA directory
cd $Nucella_netfiles/cDNA

#--------------------------------------------------------------------------------

# Concatenate reads into one fastq.gz file
#cat $Nucella_netfiles/cDNA/PROM0170_Sanford_cDNA_01102024/20240110_1712_3G_PAQ64735_daeddcc5/fastq_pass/*/*.fastq.gz \
#> $Nucella_netfiles/cDNA/Nucella.ONT.cDNA.fastq.gz

# Concatenate reads from barcode 12 (5' - GGTGCTG - CAGGTAGAAAGAAGCAGAATCGGA - TTAACCT - 3')
cat $Nucella_netfiles/cDNA/PROM0170_Sanford_cDNA_01102024/20240110_1712_3G_PAQ64735_daeddcc5/fastq_pass/barcode12/*.fastq.gz \
> $Nucella_netfiles/cDNA/Nucella.ONT.cDNA.barcode12.fastq.gz