#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=cat.reads

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --ntasks-per-node=1 # this is number of CPUs you want to use for parallel computing [also referred to as threads] 

# Reserve walltime -- hh:mm:ss
#SBATCH --time=5:00:00 #<= this may depend on your resources

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=80G #<= this may depend on your resources

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# Concatenate reads; all pass filter

# Move to the directory where you want the output files to be saved
cd /netfiles/pespenilab_share/Nucella/raw/ONT/

# Concatenate reads from each ONT cell
cat /netfiles/pespenilab_share/Nucella/raw/ONT/PROM0150_Sanford_NC3-PCR_FC2_11152023/20231115_1701_2E_PAS52172_fe64f551/fastq_pass/*.fastq.gz > FC1.ONT.nuc.fastq.gz

cat /netfiles/pespenilab_share/Nucella/raw/ONT/PROM0150_Sanford_NC3-PCR_05312023/20230531_1740_3G_PAQ05812_93506d2a/fastq_pass/*.fastq.gz > FC2.ONT.nuc.fastq.gz

cat /netfiles/pespenilab_share/Nucella/raw/ONT/PROM0172_Sanford_NC3-PCR_FC1_02282024/20240228_1655_2E_PAW14761_69937e95/fastq_pass/*.fastq.gz > FC3.ONT.nuc.fastq.gz

cat /netfiles/pespenilab_share/Nucella/raw/ONT/PROM0172_Sanford_NC3-PCR_FC2_02282024/20240228_1655_3G_PAU83574_150a145e/fastq_pass/*.fastq.gz > FC4.ONT.nuc.fastq.gz

cat /netfiles/pespenilab_share/Nucella/raw/ONT/PROM0172_Sanford_NC3-PCR_FC3_03062024/20240306_1657_1G_PAU75920_614ca090/fastq_pass/*.fastq.gz > FC5.ONT.nuc.fastq.gz

# Concatenate reads into one fastq.gz file
cat *.ONT.nuc.fastq.gz > FC_all.ONT.nuc.fastq.gz
