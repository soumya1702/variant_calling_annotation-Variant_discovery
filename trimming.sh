#!/bin/bash
#SBATCH --account=general
#SBATCH --job-name=trimming
#SBATCH --output=trimming.out
#SBATCH --error=trimming.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=72gb
#SBATCH --A students

module load trimmomatic
module load fastqc
cd /N/scratch/syennapu/GATK4

#adapter trimming
mkdir -p FASTQ/trimming/
trimmomatic PE -phred33 \
    FASTQ/S1_R1_001.fastq.gz FASTQ/S1_R2_001.fastq.gz \
    FASTQ/trimming/S1_R1_001_paired.fq.gz FASTQ/trimming/S1_R1_001_unpaired.fq.gz \
    FASTQ/trimming/S1_R2_001_paired.fq.gz FASTQ/trimming/S1_R2_001_unpaired.fq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:50
    
trimmomatic PE -phred33 \
    FASTQ/S2_R1_001.fastq.gz FASTQ/S2_R2_001.fastq.gz \
    FASTQ/trimming/S2_R1_001_paired.fq.gz FASTQ/trimming/S2_R1_001_unpaired.fq.gz \
    FASTQ/trimming/S2_R2_001_paired.fq.gz FASTQ/trimming/S2_R2_001_unpaired.fq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:50


#FASTQC
mkdir -p report/fastqc/trimfq
fastqc FASTQ/trimming/*_paired.fq.gz -o report/fastqc/trimfq/
