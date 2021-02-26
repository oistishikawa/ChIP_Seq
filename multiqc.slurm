
#!/bin/bash

#SBATCH --job-name=multiqc
#SBATCH --partition=compute
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --ntasks=1
#SBATCH --mail-user=%u@oist.jp
#SBATCH --mail-type=NONE
#SBATCH --input=none
#SBATCH --output=multiqc_%j.out
#SBATCH --error=multiqc_%j.err

module use /apps/unit/IshikawaU/.modulefiles
module load multiqc/1.9
multiqc fastq/*
