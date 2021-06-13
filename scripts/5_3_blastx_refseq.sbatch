#!/bin/bash
#SBATCH -p batch
#SBATCH -t 48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=all
#SBATCH --job-name=blastx
#SBATCH --output=stdout.blastx_refseq.%j

module load blast+/2.8.1

blastx -db refseq_protein -query Nicropus_transcriptome95.fasta -num_threads 32 -max_target_seqs 1 -out ref_seq_blast.out -outfmt '6 qseqid sseqid pident evalue bitscore stitle staxids'
