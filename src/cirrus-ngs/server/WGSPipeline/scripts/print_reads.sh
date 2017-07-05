#!/bin/bash

fileName=$1

# redirecting all output to a file
exec 1>>$2/$HOSTNAME"_printreads.o"
exec 2>>$2/$HOSTNAME"_printreads.o"

workspace=/scratch/workspace
rootdir=/shared/workspace/software

numThreads=1

bedtools=$rootdir/bedtools2/bin/bedtools
sambamba=$rootdir/sambamba/0.4.7/bin/sambamba
gatkframework=$rootdir/gatk-framework/3.2-4/bin/gatk-framework.jar
genomeFai=$rootdir/sequences/Hsapiens/ucsc.hg19.fasta.fai
genomeSeq=$rootdir/sequences/Hsapiens/ucsc.hg19.fasta

tmpDir=$workspace/tmp/$fileName

if [ ! -f $workspace/$fileName.final.bam ]; then
   java -jar -Xms166m -Xmx4g -XX:+UseSerialGC $gatkframework \
	-U LENIENT_VCF_PROCESSING \
	-T PrintReads \
	-R $genomeSeq \
	-I $workspace/$fileName.realign.bam \
	-BQSR $workspace/$fileName.dedup.grp \
	-o $workspace/$fileName.final.bam

   ### Generate reads coverage BED file ###
   $bedtools genomecov -split -ibam $workspace/$fileName.final.bam -bga -g $genomeFai -max 70001 > $workspace/$fileName.final.bed
fi

