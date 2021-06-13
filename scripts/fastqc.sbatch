#!/bin/bash
#SBATCH -p batch
#SBATCH -t 2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mail-type=all
#SBATCH --job-name=fastqc
#SBATCH --output=stdout.fastqc.%j

module load fastqc

fastqc -o fastqc_out -t 12 *.fastq.gz
