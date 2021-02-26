#!/bin/bash

#SBATCH -t 12:00:00
#SBATCH -p compute
#SBATCH --mem=100G
#SBATCH --cpus-per-task 25
#SBATCH -n 1
#SBATCH --array=0-6
#SBATCH --out=homer_%j.out

module load samtools/1.9
module load homer/v4.11

ilist_R1=(mapped/*)

input_R1=${ilist_R1[$SLURM_ARRAY_TASK_ID]}

prefix="mapped/"
output=${input_R1#"$prefix"}

makeTagDirectory TagDirectory/$output $input_R1/accepted_hits.bam
