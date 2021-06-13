#!/bin/bash
#SBATCH -p batch
#SBATCH -t 1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=all
#SBATCH --job-name=kallisto_index
#SBATCH --output=stdout.kallisto_index.%j

conda init bash
source ~/.bashrc
conda activate kallisto

kallisto index -i Nicropus_transcriptome95.fasta
