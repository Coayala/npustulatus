#!/bin/bash
#SBATCH -p batch
#SBATCH -t 120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=all
#SBATCH --job-name=hmmscan
#SBATCH --output=stdout.hmmscan.%j

module load hmmer/3.1b2

hmmscan --cpu 32 --domtblout HmmerPFAM.out Pfam-A.hmm Nicropus_transcriptome95.fasta > pfam.log
