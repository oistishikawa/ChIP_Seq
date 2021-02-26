
#!/bin/bash

#SBATCH --job-name=fastqc
#SBATCH --partition=compute
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --ntasks=1
#SBATCH --mail-user=%u@oist.jp
#SBATCH --mail-type=NONE
#SBATCH --input=none
#SBATCH --output=fastqc_%j.out
#SBATCH --error=fastqc_%j.err
#SBATCH --array=0-6

ilist=(fastq/*.fastq)
input=${ilist[$SLURM_ARRAY_TASK_ID]}

module load fastqc/0.11.5

fastqc $input
