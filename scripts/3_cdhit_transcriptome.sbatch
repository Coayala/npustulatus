#!/bin/bash
#SBATCH -p batch
#SBATCH -t 1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=all
#SBATCH --job-name=cdhit
#SBATCH --output=stdout.cdhit.%j

cd-hit-est -i Nicropus_transcriptome.fasta -o Nicropus_transcriptome95.fasta -c 0.95 -n 10 -d 0 -T 32	
