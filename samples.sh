#!/bin/bash
#SBATCH --account=general
#SBATCH --job-name=samples
#SBATCH --output=samples.out
#SBATCH --error=samples.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=72gb
#SBATCH --A students

##GATK short variants discovery pipeline

##get the dataset
mkdir -p /N/scratch/syennapu/GATK4
cd /N/scratch/syennapu/GATK4

# wget --recursive --continue --no-host-directories --no-parent --cut-dirs=1 \

tar xvf G788_spring2024/WGS_data.tar.gz

