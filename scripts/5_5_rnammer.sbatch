#!/bin/bash
#SBATCH -p batch
#SBATCH -t 120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=all
#SBATCH --job-name=rnammer
#SBATCH --output=stdout.rnammer.%j

module load rnammer/1.2

$TRINOTATE_HOME/util/rnammer_support/RnammerTranscriptome.pl --transcriptome Nicropus_transcriptome95.fasta --path_to_rnammer /opt/rnammer/1.2/prebuilt/rnammer
