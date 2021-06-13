#!/bin/bash
#SBATCH -p batch
#SBATCH -t 48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=all
#SBATCH --job-name=transdecoder
#SBATCH --output=stdout.transdecoder.%j

module load transdecoder/5.5.0

TransDecoder.LongOrfs -t Nicropus_transcriptome95.fasta
TransDecoder.Predict -t Nicropus_transcriptome95.fasta
