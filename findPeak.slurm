#!/bin/bash

#SBATCH -t 12:00:00
#SBATCH -p compute
#SBATCH --mem=100G
#SBATCH --cpus-per-task 25
#SBATCH -n 1
#SBATCH --array=0-4
#SBATCH --out=homer_%j.out

module load samtools/1.9
module load homer/v4.11

ilist_R1=(mapped/*_1*)
ilist_R2=(mapped/*_2*)

input_R1=${ilist_R1[$SLURM_ARRAY_TASK_ID]}
input_R2=${ilist_R2[$SLURM_ARRAY_TASK_ID]}

prefix="mapped/"
output=${input_R1#"$prefix"}

findPeaks TagDirectory/$sample -style factor -C 0 -o peaks/$sample.txt
