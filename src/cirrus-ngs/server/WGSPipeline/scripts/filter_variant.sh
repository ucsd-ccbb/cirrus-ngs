#!/bin/bash

workspace=/scratch/workspace
rootdir=/shared/workspace/software

export PATH=$rootdir:$PATH

threadNum=8
fileName=$1

# redirecting all output to a file
exec 1>>$2/$HOSTNAME"_filtervariants.o"
exec 2>>$2/$HOSTNAME"_filtervariants.o"

tabix=$rootdir/tabix-0.2.6/tabix
gatk=$rootdir/gatk/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar
genomeSeq=$rootdir/sequences/Hsapiens/ucsc.hg19.fasta
dbsnp=$rootdir/variation/dbsnp_138.hg19.vcf
mills=$rootdir/variation/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
omni=$rootdir/variation/1000G_omni2.5.hg19.sites.vcf
hapmap=$rootdir/variation/hapmap_3.3.hg19.sites.vcf
G1000snps=$rootdir/variation/1000G_phase1.snps.high_confidence.hg19.sites.vcf
G1000indels=$rootdir/variation/1000G_phase1.indels.hg19.sites.vcf

######################## GATK VQSR FILTRATION ###############################
# Create a Gaussian mixture model for SNPs

if [ ! -f $workspace/$fileName.vcf.gz.tbi ]; then
   tabix -p vcf $workspace/$fileName.vcf.gz
fi

if [ ! -f $workspace/$fileName.snp.tranches ]; then

java -Xmx4g -jar $gatk \
  -nt $threadNum \
  -R $genomeSeq \
  -T VariantRecalibrator \
  --maxGaussians 4 \
  -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
  -resource:omni,known=false,training=true,truth=true,prior=12.0 $omni \
  -resource:1000G,known=false,training=true,truth=false,prior=10.0 $G1000snps \
  -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
  -an QD -an MQRankSum -an ReadPosRankSum -an FS -an DP \
  -mode SNP \
  -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
  -input $workspace/$fileName.vcf.gz \
  -recalFile $workspace/$fileName.snp.recal \
  -tranchesFile $workspace/$fileName.snp.tranches \
  -rscriptFile $workspace/$fileName.snp.plots.R \
  --disable_auto_index_creation_and_locking_when_reading_rods

fi
# Create a Gaussian mixture model for INDELs

if [ ! -f $workspace/$fileName.indel.tranches ]; then

java -Xmx4g -jar $gatk \
  -nt $threadNum \
  -R $genomeSeq \
  -T VariantRecalibrator \
  --maxGaussians 4 \
  -resource:mills,known=true,training=true,truth=true,prior=12.0 $mills \
  -resource:1000G,known=false,training=true,truth=false,prior=10.0 $G1000indels \
  -an DP -an FS -an ReadPosRankSum -an MQRankSum \
  -mode INDEL \
  -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
  -input $workspace/$fileName.vcf.gz \
  -recalFile $workspace/$fileName.indel.recal \
  -tranchesFile $workspace/$fileName.indel.tranches \
  -rscriptFile $workspace/$fileName.indel.plots.R \
  --disable_auto_index_creation_and_locking_when_reading_rods

fi

if [ ! -f $workspace/$fileName.snpAr.vcf ]; then

# Apply the model for SNPs
java -Xmx4g -jar $gatk \
  -nt $threadNum \
  -R $genomeSeq \
  -T ApplyRecalibration \
  -mode SNP \
  --ts_filter_level 99.0 \
  -input $workspace/$fileName.vcf.gz \
  -recalFile $workspace/$fileName.snp.recal \
  -tranchesFile $workspace/$fileName.snp.tranches \
  -o $workspace/$fileName.snpAr.vcf

fi

if [ ! -f $workspace/$fileName.snpAr.indelAr.vcf ]; then
# Apply the model for INDELS
java -Xmx4g -jar $gatk \
  -nt $threadNum \
  -R $genomeSeq \
  -T ApplyRecalibration \
  -mode indel \
  --ts_filter_level 99.0 \
  -input $workspace/$fileName.snpAr.vcf \
  -recalFile $workspace/$fileName.indel.recal \
  -tranchesFile $workspace/$fileName.indel.tranches \
  -o $workspace/$fileName.snpAr.indelAr.vcf
fi

if [ ! -f $workspace/$fileName.vqsr.vcf ]; then
# Select the variants that passed the models
java -Xmx4g -jar $gatk \
  -nt $threadNum \
  -R $genomeSeq \
  -T SelectVariants \
  --excludeNonVariants \
  --excludeFiltered \
  --variant $workspace/$fileName.snpAr.indelAr.vcf \
  --out $workspace/$fileName.vqsr.vcf
fi
