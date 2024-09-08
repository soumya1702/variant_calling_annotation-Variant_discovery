#!/bin/bash
#SBATCH --account=general
#SBATCH --job-name=quality
#SBATCH --output=quality.out
#SBATCH --error=quality.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=72gb
#SBATCH --A students

#QC
cd /N/scratch/syennapu/GATK4

#FASTQC
module load fastqc
mkdir -p report/fastqc/rawfq
fastqc FASTQ/*.fastq.gz -o report/fastqc/rawfq/

