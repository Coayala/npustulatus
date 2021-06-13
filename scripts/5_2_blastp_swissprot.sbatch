#!/bin/bash
#SBATCH -p batch
#SBATCH -t 120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=all
#SBATCH --job-name=blastp
#SBATCH --output=stdout.blastp_swissprot.%j

module load blast+/2.8.1

blastp -query Nicropus_transcriptome95.pep -db uniprot_sprot.pep -num_threads 32 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -out blastx.out
