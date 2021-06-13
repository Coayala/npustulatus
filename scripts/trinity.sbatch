#!/bin/bash
#SBATCH -p bigmem
#SBATCH -t 120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=all
#SBATCH --job-name=trinity
#SBATCH --output=stdout.trinity.%j

module load trinity/2.8.4

Trinity --seqType fq --max_memory 750G --left ../raw_reads/FPF1/FPF1_1.fq.gz,../raw_reads/FPF2/FPF2_1.fq.gz,../raw_reads/FPF3/FPF3_1.fq.gz,../raw_reads/FPS1/FPS1_1.fq.gz,../raw_reads/FPS2/FPS2_1.fq.gz,../raw_reads/FPS3/FPS3_1.fq.gz,../raw_reads/MPF1/MPF1_1.fq.gz,../raw_reads/MPF2/MPF2_1.fq.gz,../raw_reads/MPF3/MPF3_1.fq.gz,../raw_reads/MPS1/MPS1_1.fq.gz,../raw_reads/MPS2/MPS2_1.fq.gz,../raw_reads/MPS3/MPS3_1.fq.gz  --right ../raw_reads/FPF1/FPF1_2.fq.gz,../raw_reads/FPF2/FPF2_2.fq.gz,../raw_reads/FPF3/FPF3_2.fq.gz,../raw_reads/FPS1/FPS1_2.fq.gz,../raw_reads/FPS2/FPS2_2.fq.gz,../raw_reads/FPS3/FPS3_2.fq.gz,../raw_reads/MPF1/MPF1_2.fq.gz,../raw_reads/MPF2/MPF2_2.fq.gz,../raw_reads/MPF3/MPF3_2.fq.gz,../raw_reads/MPS1/MPS1_2.fq.gz,../raw_reads/MPS2/MPS2_2.fq.gz,../raw_reads/MPS3/MPS3_2.fq.gz --CPU 32 --trimmomatic --output Nicropus_transcriptome
