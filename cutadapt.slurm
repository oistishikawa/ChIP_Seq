#!/bin/bash

#SBATCH -t 1:30:00
#SBATCH --partition=compute
#SBATCH --mem=10G
#SBATCH --cpus-per-task 32
#SBATCH --ntasks 1
#SBATCH --array=0-6
#SBATCH --output=cutadapt_%j.out
module load cutadapt/2.10

ilist=(fastq/*.fastq)

input=${ilist[$SLURM_ARRAY_TASK_ID]}

output=$(basename -s .fastq $input).fq.gz
#If the sample contains adapters,
#use the following lines to specify adapter sequence
#ADAPTER_FWD=(CTGTCTCTTATACACATCTCCGAGCCCACGAGAC)
#ADAPTER_REV=(CTGTCTCTTATACACATCTGACGCTGCCGACGA)

ADAPTER_FWD=(CTGTCTCTTATACACATCTCCGAGCCCACGAGAC)
#ADAPTER_REV=(CTGTCTCTTATACACATCTGACGCTGCCGACGA)

cutadapt -j 0 -u 40 -q 37 -g $ADAPTER_FWD --minimum-length 30 -o trimmed/$output $input
