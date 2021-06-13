#!/bin/bash
#SBATCH -p bigmem
#SBATCH -t 120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=all
#SBATCH --job-name=trinity
#SBATCH --output=stdout.trinity.%j

module load trinity/2.8.4

Trinity --seqType fq --max_memory 750G --left FPF1_1.fq.gz,FPF2_1.fq.gz,FPF3_1.fq.gz,FPS1_1.fq.gz,FPS2_1.fq.gz,FPS3_1.fq.gz,MPF1_1.fq.gz,MPF2_1.fq.gz,MPF3_1.fq.gz,MPS1_1.fq.gz,MPS2_1.fq.gz,MPS3_1.fq.gz  --right FPF1_2.fq.gz,FPF2_2.fq.gz,FPF3_2.fq.gz,FPS1_2.fq.gz,FPS2_2.fq.gz,FPS3_2.fq.gz,MPF1_2.fq.gz,MPF2_2.fq.gz,MPF3_2.fq.gz,MPS1_2.fq.gz,MPS2_2.fq.gz,MPS3_2.fq.gz --CPU 32 --trimmomatic --output Nicropus_transcriptome
