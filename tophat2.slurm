  GNU nano 2.9.8                                      tophat2.slurm                                                
#!/bin/bash

#SBATCH -t 8:00:00
#SBATCH -p compute
#SBATCH --mem=128G
#SBATCH -c 32
#SBATCH -n 2
#SBATCH --array=0-6
#SBATCH --out=tophat_%j.out

module load python/2.7.18
module load bowtie2/2.3.4.3
module load samtools/1.9
module load tophat/2.1.1

ilist=(trimmed2/*.fq.gz)

input=${ilist[$SLURM_ARRAY_TASK_ID]}

prefix="trimmed2/"
suffix=".fq.gz"
output=${input#"$prefix"}
output=${output%"$suffix"}

#echo "{$input} {$output}"

#During mapping, tophat2 generates many temp files/folders with common names
#If we run multi samples at the same time, we should seperate those temp files/folders
mkdir mapped2/$output

tophat -p32 -o mapped2/$output/ \
/flash/IshikawaU/tsunghan/mm10/Sequence/Bowtie2Index/genome \
$input
