#!/bin/bash
#SBATCH -p batch
#SBATCH -t 24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=all
#SBATCH --job-name=busco
#SBATCH --output=stdout.busco.%j

conda init bash
source ~/.bashrc

conda activate busco

busco -m transcriptome -i Nicropus_transcriptome.fasta -o Nicropus_transcriptome_busco -l arthropoda_odb10
