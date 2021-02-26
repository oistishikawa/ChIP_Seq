
#!/bin/bash
#SBATCH -p compute
#SBATCH -t 12:00:00
#SBATCH --mem=10G
#SBATCH -c 1
#SBATCH --output=sra_download_%j.out
#SBATCH --array=4795-4817

module load sra-tools/2.8.2-1

#check the SRR number which you want to download
#https://hbctraining.github.io/Accessing_public_genomic_data/lessons/downloading_from_SRA.html

# download SRR file (from SRR1124795 to SRR1124817) to specified dir and also convert it to fastq file
# if the seuquence is pair-end sequenced file, then use --split-flies

fastq-dump SRR112$SLURM_ARRAY_TASK_ID
