#!/bin/bash

workspace=/scratch/workspace
rootdir=/shared/workspace/software

export PATH=$rootdir:$PATH

fileName=$1
group=$2
# redirecting all output to a file
exec 1>>$3/$HOSTNAME"_gatkhaplotype.o"
exec 2>>$3/$HOSTNAME"_gatkhaplotype.o"

gatk=$rootdir/gatk/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar
genomeSeq=$rootdir/sequences/Hsapiens/ucsc.hg19.fasta
dbsnp138=$rootdir/variation/dbsnp_138.hg19.vcf

tmpDir=$workspace/tmp/$fileName

if [ ! -f $workspace/$fileName.raw.vcf.gz ]; then
  if [ "$group" == "NA" ]; then
     java -Xms454m -Xmx2g -XX:+UseSerialGC -Djava.io.tmpdir=$tmpDir \
	-jar $gatk \
	-T HaplotypeCaller \
	-R $genomeSeq \
	-I $workspace/$fileName.final.bam \
	-L $workspace/$fileName.final.bed \
	--annotation BaseQualityRankSumTest \
	--annotation FisherStrand \
	--annotation GCContent \
	--annotation HaplotypeScore \
	--annotation HomopolymerRun \
	--annotation MappingQualityRankSumTest \
	--annotation MappingQualityZero \
	--annotation QualByDepth \
	--annotation ReadPosRankSumTest \
	--annotation RMSMappingQuality \
	--annotation DepthPerAlleleBySample \
	--annotation Coverage \
	--annotation ClippingRankSumTest \
        --annotation DepthPerSampleHC \
        --standard_min_confidence_threshold_for_calling 30.0 \
        --standard_min_confidence_threshold_for_emitting 30.0 \
	--dbsnp $dbsnp138 \
	-o $workspace/$fileName.raw.vcf.gz \
	--disable_auto_index_creation_and_locking_when_reading_rods
   else
     ## GATK_variant_discovery_group
     java -Xms454m -Xmx7g -XX:+UseSerialGC -Djava.io.tmpdir=$tmpDir \
        -jar $gatk \
        -T HaplotypeCaller \
        -R $genomeSeq \
	-nct 2 \
        -I $workspace/$fileName.final.bam \
        -L $workspace/$fileName.final.bed \
        --annotation BaseQualityRankSumTest \
        --annotation FisherStrand \
        --annotation GCContent \
        --annotation HaplotypeScore \
        --annotation HomopolymerRun \
        --annotation MappingQualityRankSumTest \
        --annotation MappingQualityZero \
        --annotation QualByDepth \
        --annotation ReadPosRankSumTest \
        --annotation RMSMappingQuality \
        --annotation DepthPerAlleleBySample \
        --annotation Coverage \
        --annotation ClippingRankSumTest \
        --annotation DepthPerSampleHC \
        --standard_min_confidence_threshold_for_calling 30.0 \
        --standard_min_confidence_threshold_for_emitting 30.0 \
	--dbsnp $dbsnp138 \
        -o $workspace/$fileName.raw.vcf.gz \
        --disable_auto_index_creation_and_locking_when_reading_rods \
        --emitRefConfidence GVCF \
    	--variant_index_type LINEAR \
    	--variant_index_parameter 128000
   fi
fi


