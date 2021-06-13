#!/bin/bash
#SBATCH -p batch
#SBATCH -t 5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=all
#SBATCH --job-name=kallisto_quant
#SBATCH --output=stdout.kallisto_quant.%j

conda init bash
source ~/.bashrc
conda activate kallisto

# Males Fed

kallisto quant -i Nicropus_kallisto_index -t 32 -b 500 -o kallisto_Nicropus-counts/MPF1/ -t 32 MPF1_1.fq.gz MPF1_2.fq.gz
kallisto quant -i Nicropus_kallisto_index -t 32 -b 500 -o kallisto_Nicropus-counts/MPF2/ -t 32 MPF2_1.fq.gz MPF2_2.fq.gz
kallisto quant -i Nicropus_kallisto_index -t 32 -b 500 -o kallisto_Nicropus-counts/MPF3/ -t 32 MPF3_1.fq.gz MPF3_2.fq.gz

# Males starved

kallisto quant -i Nicropus_kallisto_index -t 32 -b 500 -o kallisto_Nicropus-counts/MPS1/ -t 32 MPS1_1.fq.gz MPS1_2.fq.gz
kallisto quant -i Nicropus_kallisto_index -t 32 -b 500 -o kallisto_Nicropus-counts/MPS2/ -t 32 MPS2_1.fq.gz MPS2_2.fq.gz
kallisto quant -i Nicropus_kallisto_index -t 32 -b 500 -o kallisto_Nicropus-counts/MPS3/ -t 32 MPS3_1.fq.gz MPS3_2.fq.gz

# Females Fed

kallisto quant -i Nicropus_kallisto_index -t 32 -b 500 -o kallisto_Nicropus-counts/FPF1/ -t 32 FPF1_1.fq.gz FPF1_2.fq.gz
kallisto quant -i Nicropus_kallisto_index -t 32 -b 500 -o kallisto_Nicropus-counts/FPF2/ -t 32 FPF2_1.fq.gz FPF2_2.fq.gz
kallisto quant -i Nicropus_kallisto_index -t 32 -b 500 -o kallisto_Nicropus-counts/FPF3/ -t 32 FPF3_1.fq.gz FPF3_2.fq.gz

# Females starved

kallisto quant -i Nicropus_kallisto_index -t 32 -b 500 -o kallisto_Nicropus-counts/FPS1/ -t 32 FPS1_1.fq.gz FPS1_2.fq.gz
kallisto quant -i Nicropus_kallisto_index -t 32 -b 500 -o kallisto_Nicropus-counts/FPS2/ -t 32 FPS2_1.fq.gz FPS2_2.fq.gz
kallisto quant -i Nicropus_kallisto_index -t 32 -b 500 -o kallisto_Nicropus-counts/FPS3/ -t 32 FPS3_1.fq.gz FPS3_2.fq.gz