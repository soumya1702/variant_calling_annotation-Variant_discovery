#!/bin/bash
#SBATCH --account=general
#SBATCH --job-name=alignment
#SBATCH --output=alignment.out
#SBATCH --error=alignment.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=72gb
#SBATCH --A students


#alignment for 2 samples s1 ans s2
module load bwa
module load gatk
module load samtools
module load picard
picard=/N/soft/rhel8/picard/3.0.0/build/libs/picard.jar

#alignment
GATK_bundle=Reference/broad.hg38.v0/
reference=$GATK_bundle/Homo_sapiens_assembly38.fasta
NT=8
SID=S1
# SID=S2
read_group="@RG\tID:${SID}\tSM:${SID}\tPL:ILLUMINA"

mkdir -p report/alignment/${SID}/
bwa mem \
    -t $NT \
    -R ${read_group} \
    ${reference} \
    FASTQ/trimming/${SID}_R1_001_paired.fq.gz \
    FASTQ/trimming/${SID}_R2_001_paired.fq.gz \
    > report/alignment/${SID}/${SID}.sam

#sorting and indexing
samtools view -bS report/alignment/${SID}/${SID}.sam > report/alignment/${SID}/${SID}.bam 
samtools sort -@ $NT -o report/alignment/${SID}/${SID}.sorted.bam report/alignment/${SID}/${SID}.bam 
samtools index report/alignment/${SID}/${SID}.sorted.bam

#mark duplicates
module load picard
picard=/N/soft/rhel8/picard/3.0.0/build/libs/picard.jar
java -jar $picard MarkDuplicates \
    I=report/alignment/${SID}/${SID}.sorted.bam \
    O=report/alignment/${SID}/${SID}.sorted.mkdp.bam \
    M=report/alignment/${SID}/${SID}.sorted.mkdp_metrics.txt

samtools index report/alignment/${SID}/${SID}.sorted.mkdp.bam

#================================================================================
#BQSR
SID=S1
# SID=S2

gatk BaseRecalibrator \
    -R ${reference} \
    -I report/alignment/${SID}/${SID}.sorted.mkdp.bam \
    --known-sites $GATK_bundle/Homo_sapiens_assembly38.dbsnp138.vcf \
    --known-sites $GATK_bundle/Homo_sapiens_assembly38.known_indels.vcf.gz \
    -O report/alignment/${SID}/${SID}.recal_data.table
 
gatk ApplyBQSR \
   -R ${reference} \
   -I report/alignment/${SID}/${SID}.sorted.mkdp.bam \
   --bqsr-recal-file report/alignment/${SID}/${SID}.recal_data.table \
   -O report/alignment/${SID}/${SID}.sorted.mkdp.bqsr.bam
############################################################################################
##########################################################################################
#===========================================================================================
#alignment
GATK_bundle=Reference/broad.hg38.v0/
reference=$GATK_bundle/Homo_sapiens_assembly38.fasta
NT=8
#SID=S1
SID=S2
read_group="@RG\tID:${SID}\tSM:${SID}\tPL:ILLUMINA"

mkdir -p report/alignment/${SID}/
bwa mem \
    -t $NT \
    -R ${read_group} \
    ${reference} \
    FASTQ/trimming/${SID}_R1_001_paired.fq.gz \
    FASTQ/trimming/${SID}_R2_001_paired.fq.gz \
    > report/alignment/${SID}/${SID}.sam

#sorting and indexing
samtools view -bS report/alignment/${SID}/${SID}.sam > report/alignment/${SID}/${SID}.bam 
samtools sort -@ $NT -o report/alignment/${SID}/${SID}.sorted.bam report/alignment/${SID}/${SID}.bam 
samtools index report/alignment/${SID}/${SID}.sorted.bam

#mark duplicates
module load picard
picard=/N/soft/rhel8/picard/3.0.0/build/libs/picard.jar
java -jar $picard MarkDuplicates \
    I=report/alignment/${SID}/${SID}.sorted.bam \
    O=report/alignment/${SID}/${SID}.sorted.mkdp.bam \
    M=report/alignment/${SID}/${SID}.sorted.mkdp_metrics.txt

samtools index report/alignment/${SID}/${SID}.sorted.mkdp.bam

#================================================================================
#BQSR
#SID=S1
SID=S2

gatk BaseRecalibrator \
    -R ${reference} \
    -I report/alignment/${SID}/${SID}.sorted.mkdp.bam \
    --known-sites $GATK_bundle/Homo_sapiens_assembly38.dbsnp138.vcf \
    --known-sites $GATK_bundle/Homo_sapiens_assembly38.known_indels.vcf.gz \
    -O report/alignment/${SID}/${SID}.recal_data.table
 
gatk ApplyBQSR \
   -R ${reference} \
   -I report/alignment/${SID}/${SID}.sorted.mkdp.bam \
   --bqsr-recal-file report/alignment/${SID}/${SID}.recal_data.table \
   -O report/alignment/${SID}/${SID}.sorted.mkdp.bqsr.bam

