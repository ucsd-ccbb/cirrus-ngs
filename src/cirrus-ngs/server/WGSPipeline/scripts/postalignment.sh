#!/bin/bash

fileName=$1

# redirecting all output to a file
exec 1>>$2/$HOSTNAME"_postalignment.o"
exec 2>>$2/$HOSTNAME"_postalignment.o"

workspace=/scratch/workspace
rootdir=/shared/workspace/software

numThreads=1

samtools=$rootdir/samtools/samtools-1.1/samtools
bedtools=$rootdir/bedtools2/bin/bedtools
sambamba=$rootdir/sambamba/0.4.7/bin/sambamba
gatk=$rootdir/gatk/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar
picard=$rootdir/picard-1.96
genomeFai=$rootdir/sequences/Hsapiens/ucsc.hg19.fasta.fai
genomeSeq=$rootdir/sequences/Hsapiens/ucsc.hg19.fasta
dbsnp=$rootdir/variation/dbsnp_138.hg19.vcf
mills=$rootdir/variation/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
G1000=$rootdir/variation/1000G_phase1.snps.high_confidence.hg19.sites.vcf

tmpDir=$workspace/tmp/$fileName

if [ ! -f $workspace/$fileName.dedup.bam ]; then
   java -jar -Xms250m -Xmx2g $picard/MarkDuplicates.jar \
	INPUT=$workspace/$fileName.sort.split.bam \
	OUTPUT=$workspace/$fileName.dedup.bam \
	METRICS_FILE=$workspace/$fileName.metrics.txt \
	AS=true VALIDATION_STRINGENCY=LENIENT

   $sambamba index -t $numThreads $workspace/$fileName.dedup.bam $workspace/$fileName.dedup.bam.bai
fi

### Generate reads coverage BED file ###
if [ ! -f $workspace/$fileName.dedup.bed ]; then
   $bedtools genomecov -split -ibam $workspace/$fileName.dedup.bam -bga -g $genomeFai -max 70001 > $workspace/$fileName.dedup.bed
fi

### RealignerTargetCreator ###
if [ ! -f $workspace/$fileName.realign.intervals ]; then
   java -Xms250m -Xmx2g -XX:+UseSerialGC -Djava.io.tmpdir=$tmpDir \
	-jar $gatk \
	-T RealignerTargetCreator \
	-I $workspace/$fileName.dedup.bam \
	-R $genomeSeq \
	-o $workspace/$fileName.realign.intervals \
	-l INFO --interval_set_rule INTERSECTION \
	--known $mills \
	--known $dbsnp \
        --read_filter BadCigar \
        --read_filter NotPrimaryAlignment \
	-U LENIENT_VCF_PROCESSING \
	--disable_auto_index_creation_and_locking_when_reading_rods
fi

### IndelRealigner ###
if [ ! -f $workspace/$fileName.realign.bam ]; then
   java -Xms454m -Xmx2g -XX:+UseSerialGC -Djava.io.tmpdir=$tmpDir \
	-jar $gatk \
	-T IndelRealigner \
	-I $workspace/$fileName.dedup.bam \
	-R $genomeSeq \
	-targetIntervals $workspace/$fileName.realign.intervals \
	-U LENIENT_VCF_PROCESSING \
	-known $mills \
	-known $dbsnp \
	--disable_auto_index_creation_and_locking_when_reading_rods \
	--read_filter BadCigar \
	--read_filter NotPrimaryAlignment \
	-o $workspace/$fileName.realign.bam
fi

### BaseRecalibrator ###
if [ ! -f $workspace/$fileName.dedup.grp ]; then
   java -Xms454m -Xmx2g -XX:+UseSerialGC -Djava.io.tmpdir=$tmpDir \
	-jar $gatk \
	-T BaseRecalibrator \
	-o $workspace/$fileName.dedup.grp \
	-I $workspace/$fileName.dedup.bam \
	-L $workspace/$fileName.dedup.bed \
	-R $genomeSeq \
	--knownSites $mills \
	--knownSites $dbsnp \
	--knownSites $G1000 \
	--interval_set_rule INTERSECTION \
	-U LENIENT_VCF_PROCESSING \
	--read_filter BadCigar \
        --read_filter NotPrimaryAlignment \
	--disable_auto_index_creation_and_locking_when_reading_rods
fi

