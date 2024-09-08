#!/bin/bash
#SBATCH --account=general
#SBATCH --job-name=haplotypecaller
#SBATCH --output=haplotypecaller.out
#SBATCH --error=haplotypecaller.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=72gb
#SBATCH --A students


module load gatk
module load picard
module load samtools
picard=/N/soft/rhel8/picard/3.0.0/build/libs/picard.jar

##variant calling for large scale study

##Run alignment and BQSR for all samples
##Run the code above before "HaplotypeCaller" with SID=S2

##Run HaplotypeCaller on each sample
for SID in S1 S2
do
    if [[ ! -f report/gatk/${SID}/${SID}.g.vcf.gz ]];then 
        gatk HaplotypeCaller \
          -R ${reference} \
          -I report/alignment/${SID}/${SID}.sorted.mkdp.bqsr.bam \
          -ERC GVCF \
          -O report/gatk/${SID}/${SID}.g.vcf.gz
    fi
done

##run by chromosome separately. this dataset has only chr11 reads
test -e report/gatk/GenomicsDB && rm -rf report/gatk/GenomicsDB #delete the folder if it exists
gatk GenomicsDBImport \
    -V report/gatk/S1/S1.g.vcf.gz \
    -V report/gatk/S2/S2.g.vcf.gz \
    -L chr11 \
    --genomicsdb-workspace-path report/gatk/GenomicsDB

gatk GenotypeGVCFs \
    -R ${reference} \
    -V gendb://report/gatk/GenomicsDB \
    -O report/gatk/joint_genotyped.vcf.gz

##Run GATK VariantRecalibrator for SNPs
input_vcf=report/gatk/joint_genotyped.vcf.gz
mkdir -p report/gatk/VQSR

gatk VariantRecalibrator \
    -R $reference \
    -V $input_vcf \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $GATK_bundle/hapmap_3.3.hg38.vcf.gz \
    -resource:omni,known=false,training=true,truth=false,prior=12.0 $GATK_bundle/1000G_omni2.5.hg38.vcf.gz \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 $GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $GATK_bundle/Homo_sapiens_assembly38.dbsnp138.vcf \
    -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum \
    -mode SNP \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
    --tranches-file report/gatk/VQSR/snp.tranches \
    --rscript-file report/gatk/VQSR/output_SNP.plots.R \
    -O report/gatk/VQSR/snp.recal

##Apply VQSR to the SNPs in the VCF
gatk ApplyVQSR \
    -R $reference \
    -V $input_vcf \
    --recal-file report/gatk/VQSR/snp.recal \
    --tranches-file report/gatk/VQSR/snp.tranches \
    --truth-sensitivity-filter-level 99.0 \
    -mode SNP \
    -O report/gatk/VQSR/snp.output.vcf

##Run GATK VariantRecalibrator for INDELs
gatk VariantRecalibrator \
    -R $reference \
    -V report/gatk/VQSR/snp.output.vcf \
    -resource:mills,known=false,training=true,truth=true,prior=12.0 $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $GATK_bundle/Homo_sapiens_assembly38.known_indels.vcf.gz \
    -an DP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum \
    -mode INDEL \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
    --tranches-file report/gatk/VQSR/indel.tranches \
    --rscript-file report/gatk/VQSR/output_indel.plots.R \
    -O report/gatk/VQSR/indel.recal

##Apply VQSR to the INDELs in the VCF
gatk ApplyVQSR \
    -R $reference \
    -V report/gatk/VQSR/snp.output.vcf \
    --recal-file report/gatk/VQSR/indel.recal \
    --tranches-file report/gatk/VQSR/indel.tranches \
    --truth-sensitivity-filter-level 99.0 \
    -mode INDEL \
    -O report/gatk/final.output.vcf


