#!/bin/bash

fileName=$1
chromosome=$2

workspace=/scratch/workspace/$fileName
#workspace=/shared/shuling/data
#gatk=/shared/workspace/software/gatk/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar
gatk=/shared/shuling/GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar

aws s3 cp s3://analysis-data-by-ccbb/20170320_Shuling_Califano_dnaseq/analysis_results/$fileName/$fileName.final.$chromosome.bam $workspace/
aws s3 cp s3://analysis-data-by-ccbb/20170320_Shuling_Califano_dnaseq/analysis_results/$fileName/$fileName.final.$chromosome.bai $workspace/

mkdir -p $workspace/temp

# redirecting all output to a file
exec 1>>$workspace/"variant_call.o"
exec 2>>$workspace/"variant_call.e"

/shared/workspace/software/bedtools2/bin/bedtools genomecov -split -ibam $workspace/$fileName.final.$chromosome.bam -bga -g /shared/shuling/resource/bwa/human_g1k_v37.fasta.fai -max 70001 > $workspace/$fileName.final.$chromosome.bed

java -Xms454m -Xmx8g -XX:+UseSerialGC -Djava.io.tmpdir=$workspace/temp \
  -jar $gatk \
  -T HaplotypeCaller \
  -R /shared/shuling/resource/bwa/human_g1k_v37.fasta \
  -I $workspace/$fileName.final.$chromosome.bam \
  -L $workspace/$fileName.final.$chromosome.bed \
  --out $workspace/$fileName.raw.$chromosome.vcf.gz \
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
  --standard_min_confidence_threshold_for_calling 30.0 \
  --standard_min_confidence_threshold_for_emitting 30.0 \
  --dbsnp /shared/shuling/resource/dbsnp_132_b37.leftAligned.vcf


#  --annotation DepthPerSampleHC \

aws s3 cp $workspace/$fileName.raw.$chromosome.vcf.gz s3://analysis-data-by-ccbb/20170320_Shuling_Califano_dnaseq/analysis_results/variants/$fileName/

