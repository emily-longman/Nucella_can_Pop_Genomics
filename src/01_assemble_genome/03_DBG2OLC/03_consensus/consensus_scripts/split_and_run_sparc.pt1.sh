#!/usr/bin/env bash

###
# USAGE: ./split_and_run_sparc.sh [BACKBONE_FASTA] [CONSENSUS_FASTA] [READS_FASTA] [OUTPUT_DIR] [ITERATIONS] ###

backbone_fasta=$1 #backbone broken into chunks
consensus_fasta=$2 #text file with contig names in chunks
reads_fasta=$3 #contigs and raw ONT reads catted together
split_dir=$4 #file output directory 
iterations=$5 #set to 2
ncpus=$6 #set to 32
nshift=$7 #shiftn

echo $backbone_fasta
echo $consensus_fasta
echo $reads_fasta
echo $split_dir
echo $iterations
echo $ncpus
echo $nshift


### need scripts
##split_reads_by_backbone.py

### Make some aliaces
#alias blasr="/gpfs1/home/j/c/jcnunez/software/blasrmc/alignment/bin/blasrmc" #now loaded from the VACC
alias Sparc="/gpfs1/home/e/l/elongman/software/Sparc"

#clean the directory first  --- no longer required
#find ${split_dir} -name "backbone-*" -delete

/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/00_Genome_short_read/12_consensus_scripts_extra/split_reads_by_backbone.py -b ${backbone_fasta} -o ${split_dir} -r ${reads_fasta} -c ${consensus_fasta} -s ${nshift}

#===#for file in $(find ${split_dir} -name "*.reads.fasta"); do
#===#chunk=`basename $file .reads.fasta`
#===#
#===#cmd=""
#===#for iter in `seq 1 ${iterations}`; do
#===## create this to shorten the subsequent commands to more readable length
#===## I think I got them all substituted correctly.
#===#d="${split_dir}/${chunk}"
#===#cmd="blasr -nproc $ncpus ${d}.reads.fasta ${d}.fasta -bestn 1 -m 5 -minMatch 19 -out ${d}.mapped.m5"
#===#echo $cmd ; eval $cmd
#===#cmd="Sparc m ${d}.mapped.m5 b ${d}.fasta k 1 c 2 g 1 HQ_Prefix Contig boost 5 t 0.2 o ${d}"
#===#echo $cmd ; eval $cmd        
#===#if [[ ${iter} -lt ${iterations} ]]
#===#then
#===## rename
#===#cmd="mv ${d}.consensus.fasta ${d}.fasta"
#===#echo $cmd ; eval $cmd
#===#fi
#===#done
#===#
#===#echo $cmd
#===#eval $cmd
#===#
#===#
#===##to save space
#===#cmd="rm ${split_dir}/${chunk}.mapped.m5"
#===#echo $cmd
#===#eval $cmd
#===#cmd="rm ${split_dir}/${chunk}.reads.fasta"
#===#echo $cmd
#===#eval $cmd
#===#
#===#done
#===#
#===#for confile in $(find ${split_dir} -name "*.consensus.fasta"); do
#===#cmd="cat ${confile};"
#===#eval $cmd
#===#done > /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly/consensus/final_assembly.fasta