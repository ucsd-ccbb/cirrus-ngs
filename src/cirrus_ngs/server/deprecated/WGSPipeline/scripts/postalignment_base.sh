#!/bin/bash

fileName=$1

for chromosome in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT
do

workspace=/scratch/workspace/$fileName

mkdir $workspace/logs/

aws s3 cp s3://analysis-data-by-ccbb/20170320_Shuling_Califano_dnaseq/analysis_results/$fileName/$fileName.$chromosome.dedup.bam $workspace/
aws s3 cp s3://analysis-data-by-ccbb/20170320_Shuling_Califano_dnaseq/analysis_results/$fileName/$fileName.$chromosome.dedup.bam.bai $workspace/
aws s3 cp s3://analysis-data-by-ccbb/20170320_Shuling_Califano_dnaseq/analysis_results/$fileName/$fileName.$chromosome.dedup.bed $workspace/
aws s3 cp s3://analysis-data-by-ccbb/20170320_Shuling_Califano_dnaseq/analysis_results/$fileName/$fileName.realign.$chromosome.bam $workspace/
aws s3 cp s3://analysis-data-by-ccbb/20170320_Shuling_Califano_dnaseq/analysis_results/$fileName/$fileName.realign.$chromosome.bai $workspace/

if [ ! -f $workspace/$fileName.$chromosome.grp ]; then
echo "$(date): running BaseRecalibrator..."
## BaseRecalibrator
java -Djava.io.tmpdir=$workspace/temp -Xmx8g -jar /shared/shuling/GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar -T BaseRecalibrator -R /shared/shuling/resource/bwa/human_g1k_v37.fasta -I $workspace/$fileName.$chromosome.dedup.bam -knownSites /shared/shuling/resource/dbsnp_138.b37.vcf -knownSites /shared/shuling/resource/hapmap_3.3.b37.vcf -o $workspace/$fileName.$chromosome.grp -dcov 2 -L $workspace/$fileName.$chromosome.dedup.bed --allow_potentially_misencoded_quality_scores
fi

if [ ! -f $workspace/$fileName.final.$chromosome.bam ]; then
echo "$(date): running PrintReads..."
## PrintReads
java -Djava.io.tmpdir=$workspace/temp -Xmx8g -jar /shared/shuling/GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar -T PrintReads -R /shared/shuling/resource/bwa/human_g1k_v37.fasta -I $workspace/$fileName.realign.$chromosome.bam -BQSR $workspace/$fileName.$chromosome.grp -o $workspace/$fileName.final.$chromosome.bam -rf BadCigar -L $chromosome --no_pg_tag --allow_potentially_misencoded_quality_scores
fi

aws s3 cp $workspace/$fileName.$chromosome.grp s3://analysis-data-by-ccbb/20170320_Shuling_Califano_dnaseq/analysis_results/$fileName/
aws s3 cp $workspace/$fileName.final.$chromosome.bam s3://analysis-data-by-ccbb/20170320_Shuling_Califano_dnaseq/analysis_results/$fileName/
aws s3 cp $workspace/$fileName.final.$chromosome.bai s3://analysis-data-by-ccbb/20170320_Shuling_Califano_dnaseq/analysis_results/$fileName/

done
