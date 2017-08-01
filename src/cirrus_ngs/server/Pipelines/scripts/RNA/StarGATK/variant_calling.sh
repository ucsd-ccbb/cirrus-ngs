#!/bin/bash

rootdir=/shared/workspace/software
picard=$rootdir/picard-1.96
gatk=$rootdir/gatk/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar
sambamba=$rootdir/sambamba/0.4.7/bin/sambamba
bgzip=$rootdir/tabix-0.2.6/bgzip
tabix=$rootdir/tabix-0.2.6/tabix
refFasta=$rootdir/sequences/Hsapiens/ucsc.hg19.fasta
mills=$rootdir/variation/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
g1000=$rootdir/variation/1000G_phase1.snps.high_confidence.hg19.sites.vcf
dbsnp=$rootdir/variation/dbsnp_138.hg19.vcf

# create output dir
sample=$1
project_name=$2
workflow=$3

workspace=/shared/workspace/data_archive/RNASeq/$project_name/$workflow
inDir=$workspace/$sample
outDir=$workspace/$sample/variants
tempDir=$workspace/tmp

mkdir -p $outDir/$tempDir

java -Djava.io.tmpdir=$tempDir -jar $picard/AddOrReplaceReadGroups.jar \
        I=$inDir/Aligned.out.bam \
        O=$outDir/Aligned.sortedByCoord.out.grp.bam \
	SO=coordinate \
	RGID=ILLUMINA RGLB=LaneX RGPL=illumina RGPU=NONE RGSM=HLI_tumor

# MARK DUPLICATES
java -jar $picard/MarkDuplicates.jar \
	I=$outDir/Aligned.sortedByCoord.out.grp.bam \
	O=$outDir/Aligned.sortedByCoord.out.dedupped.bam \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=SILENT \
	M=$outDir/output.metrics

# SPLIT'N'TRIM + MAPPING Q
java -Djava.io.tmpdir=$tempDir -Xmx15g -jar $gatk \
	-T SplitNCigarReads \
	-R $refFasta \
	-I $outDir/Aligned.sortedByCoord.out.dedupped.bam \
	-o $outDir/Aligned.sortedByCoord.out.dedupped.split.bam \
	-rf ReassignOneMappingQuality \
	-RMQF 255 \
	-RMQT 60 \
	-U ALLOW_N_CIGAR_READS

# INDEL REALIGNMENT
java -Djava.io.tmpdir=$tempDir -Xmx15g -jar $gatk \
	-T RealignerTargetCreator \
	-R $refFasta \
	-I $outDir/Aligned.sortedByCoord.out.dedupped.split.bam \
	-o $outDir/output.intervals \
	-known $mills \
	-known $g1000 \
	-nt 8

java -Djava.io.tmpdir=$tempDir -Xmx15g -jar $gatk \
	-I $outDir/Aligned.sortedByCoord.out.dedupped.split.bam \
	-R $refFasta \
	-T IndelRealigner \
	-targetIntervals $outDir/output.intervals \
	-o $outDir/Aligned.sortedByCoord.out.dedupped.split.realigned.bam \
	-known $mills \
	-known $g1000 \
	--consensusDeterminationModel KNOWNS_ONLY \
	-LOD 0.4

# BQSR
java -Djava.io.tmpdir=$tempDir -Xmx15g -jar $gatk \
    -T BaseRecalibrator \
    -R $refFasta \
    -I $outDir/Aligned.sortedByCoord.out.dedupped.split.realigned.bam \
    -knownSites $dbsnp \
    -knownSites $mills \
    -knownSites $g1000 \
    -o $outDir/recal_data.table \
    -nct 10

java -Djava.io.tmpdir=$tempDir -Xmx15g -jar $gatk \
    -T PrintReads \
    -R $refFasta \
    -I $outDir/Aligned.sortedByCoord.out.dedupped.split.realigned.bam \
    -BQSR $outDir/recal_data.table \
    -o $outDir/Aligned.sortedByCoord.out.dedupped.split.realigned.recal.bam

# VARIANT CALLING
java -Djava.io.tmpdir=$tempDir -Xmx15g -jar $gatk \
	-T HaplotypeCaller \
	-R $refFasta \
	-I $outDir/Aligned.sortedByCoord.out.dedupped.split.realigned.recal.bam \
	-dontUseSoftClippedBases \
	-stand_call_conf 20.0 \
	-stand_emit_conf 20.0 \
	-o $outDir/${sample}_raw.vcf

# VARIANT FILTERING
java -Djava.io.tmpdir=$tempDir -Xmx15g -jar $gatk \
	-T VariantFiltration \
	-R $refFasta \
	-V $outDir/${sample}_raw.vcf \
	-window 35 \
	-cluster 3 \
	-filterName FS \
	-filter "FS > 30.0" \
	-filterName QD \
	-filter "QD < 2.0" \
	-o $outDir/${sample}_filt.vcf

# zip and index
$bgzip $outDir/${sample}_filt.vcf -c > $outDir/${sample}_filt.vcf.gz
$tabix -p vcf $outDir/${sample}_filt.vcf.gz


