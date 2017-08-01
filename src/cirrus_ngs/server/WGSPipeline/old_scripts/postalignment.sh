#!/bin/bash

fileName=$1
chromosome=$2
workspace=/scratch/workspace/$fileName

mkdir $workspace/logs/

# redirecting all output to a file
#exec 1>>$workspace/logs/"postalignment.o"
#exec 2>>$workspace/logs/"postalignment.e"

if [ ! -f $workspace/$fileName.$chromosome.dedup.bam ]; then
echo "$(date): running picard..."
java -jar -Djava.io.tmpdir=$workspace/temp -Xms250m -Xmx20g /shared/workspace/software/picard-1.96/MarkDuplicates.jar INPUT=$workspace/$fileName.$chromosome.bam OUTPUT=$workspace/$fileName.$chromosome.dedup.bam METRICS_FILE=$workspace/$fileName.metrics.txt AS=true VALIDATION_STRINGENCY=LENIENT
fi

if [ ! -f $workspace/$fileName.$chromosome.dedup.bam.bai ]; then
echo "$(date): running sambamba..."
## Sambamba
/shared/workspace/software/sambamba/0.4.7/bin/sambamba index -t 1 $workspace/$fileName.$chromosome.dedup.bam $workspace/$fileName.$chromosome.dedup.bam.bai
fi

if [ ! -f $workspace/$fileName.$chromosome.dedup.bed ]; then
echo "$(date): running bedtools..."
## Bedtools
/shared/workspace/software/bedtools2/bin/bedtools genomecov -split -ibam $workspace/$fileName.$chromosome.dedup.bam -bga -g /shared/shuling/resource/bwa/human_g1k_v37.fasta.fai -max 70001 > $workspace/$fileName.$chromosome.dedup.bed
fi

if [ ! -f $workspace/intervals.$chromosome.interval_list ]; then
echo "$(date): running RealignerTargetCreator..."
## RealignerTargetCreator
java -Djava.io.tmpdir=$workspace/temp -Xmx4g -jar /shared/shuling/GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /shared/shuling/resource/bwa/human_g1k_v37.fasta -I $workspace/$fileName.$chromosome.dedup.bam -known /shared/shuling/resource/1000G_phase1.indels.b37.vcf -known /shared/shuling/resource/Mills_and_1000G_gold_standard.indels.b37.vcf -o $workspace/intervals.$chromosome.interval_list -L $chromosome --allow_potentially_misencoded_quality_scores
fi

if [ ! -f $workspace/$fileName.realign.$chromosome.bam ]; then
echo "$(date): running IndelRealigner..."
java -Djava.io.tmpdir=$workspace/temp -Xmx4g -jar /shared/shuling/GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar -T IndelRealigner -R /shared/shuling/resource/bwa/human_g1k_v37.fasta -I $workspace/$fileName.$chromosome.dedup.bam -known /shared/shuling/resource/1000G_phase1.indels.b37.vcf -known /shared/shuling/resource/Mills_and_1000G_gold_standard.indels.b37.vcf -targetIntervals $workspace/intervals.$chromosome.interval_list -o $workspace/$fileName.realign.$chromosome.bam -L $chromosome --allow_potentially_misencoded_quality_scores
fi

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

