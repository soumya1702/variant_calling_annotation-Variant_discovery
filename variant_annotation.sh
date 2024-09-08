#!/bin/bash
#SBATCH --account=general
#SBATCH --job-name=variant_annotation
#SBATCH --output=variant_annotation.out
#SBATCH --error=variant_annotation.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=72gb
#SBATCH --A students

module load gatk
module load picard
module load samtools

##variant annotation using ANNOVAR

##for one-sample analysis
#gunzip report/gatk/S1/S1.hard.vcf.gz -c > report/gatk/S1/S1.hard.vcf
#input_vcf=report/gatk/S1/S1.hard.vcf

##for multi-sample analysis
input_vcf=report/gatk/final.output.vcf

mkdir report/ANNOVAR

table_annovar.pl \
  $input_vcf \
  hg38/ \
  -buildver hg38 \
  -out report/ANNOVAR/WGS \
  -remove \
  -protocol refGene,cytoBand,exac03,avsnp150,dbnsfp41a \
  -operation g,r,f,f,f \
  -nastring . \
  -vcfinput
