#!/bin/bash

# rootdir=/shared/workspace/software
# picard=$rootdir/picard-1.96        # TODO: update
# gatk=$rootdir/gatk/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar
# sambamba=$rootdir/sambamba/0.4.7/bin/sambamba
# bgzip=$rootdir/tabix-0.2.6/bgzip
# tabix=$rootdir/tabix-0.2.6/tabix
# genome_fasta=$rootdir/sequences/Hsapiens/ucsc.hg19.fasta
# mills=$rootdir/variation/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
# G1000=$rootdir/variation/1000G_phase1.snps.high_confidence.hg19.sites.vcf
# dbsnp=$rootdir/variation/dbsnp_138.hg19.vcf

project_name=$1
file_suffix=$2  #extension of input file, does not include .gz if present in input
root_dir=$3
fastq_end1=$4
fastq_end2=$5
input_address=$6    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$7   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$8
is_zipped=$9    #either "True" or "False", indicates whether input is gzipped
num_threads=${10}       #number of threads

#logging
log_dir=$log_dir/$fastq_end1
mkdir -p $log_dir
log_file=$log_dir/'gatk_vc.log'
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name
mkdir -p $workspace

inDir=$workspace/$fastq_end1
outDir=$workspace/$fastq_end1/variants
tempDir=$workspace/tmp

mkdir -p $outDir/$tempDir

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME $status_file

##DOWNLOAD##
if [ ! -f $inDir/$fastq_end1$file_suffix ]
then
    aws s3 cp $input_address/$fastq_end1$file_suffix $inDir/
fi
##END_DOWNLOAD##

# START THE PIPELINE

check_exit_status "$java -Djava.io.tmpdir=$tempDir -jar $picard/AddOrReplaceReadGroups.jar \
    I=$inDir/$fastq_end1$file_suffix \
    O=$outDir/Aligned.sortedByCoord.out.grp.bam \
	SO=coordinate \
	RGID=ILLUMINA RGLB=LaneX RGPL=illumina RGPU=NONE RGSM=HLI_tumor" $JOB_NAME $status_file

# MARK DUPLICATES
check_exit_status "$java -jar $picard/MarkDuplicates.jar \
	I=$outDir/Aligned.sortedByCoord.out.grp.bam \
	O=$outDir/Aligned.sortedByCoord.out.dedupped.bam \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=SILENT \
	M=$outDir/output.metrics" $JOB_NAME $status_file

# SPLIT'N'TRIM + MAPPING Q
check_exit_status "$java -Djava.io.tmpdir=$tempDir -Xmx15g -jar $gatk \
	-T SplitNCigarReads \
	-R $genome_fasta \
	-I $outDir/Aligned.sortedByCoord.out.dedupped.bam \
	-o $outDir/Aligned.sortedByCoord.out.dedupped.split.bam \
	-rf ReassignOneMappingQuality \
	-RMQF 255 \
	-RMQT 60 \
	-U ALLOW_N_CIGAR_READS" $JOB_NAME $status_file

# INDEL REALIGNMENT
check_exit_status "$java -Djava.io.tmpdir=$tempDir -Xmx15g -jar $gatk \
	-T RealignerTargetCreator \
	-R $genome_fasta \
	-I $outDir/Aligned.sortedByCoord.out.dedupped.split.bam \
	-o $outDir/output.intervals \
	-known $mills \
	-known $G1000 \
	-nt 8" $JOB_NAME $status_file

check_exit_status "$java -Djava.io.tmpdir=$tempDir -Xmx15g -jar $gatk \
	-I $outDir/Aligned.sortedByCoord.out.dedupped.split.bam \
	-R $genome_fasta \
	-T IndelRealigner \
	-targetIntervals $outDir/output.intervals \
	-o $outDir/Aligned.sortedByCoord.out.dedupped.split.realigned.bam \
	-known $mills \
	-known $G1000 \
	--consensusDeterminationModel KNOWNS_ONLY \
	-LOD 0.4" $JOB_NAME $status_file

# BQSR
check_exit_status "$java -Djava.io.tmpdir=$tempDir -Xmx15g -jar $gatk \
    -T BaseRecalibrator \
    -R $genome_fasta \
    -I $outDir/Aligned.sortedByCoord.out.dedupped.split.realigned.bam \
    -knownSites $dbsnp \
    -knownSites $mills \
    -knownSites $G1000 \
    -o $outDir/recal_data.table \
    -nct $num_threads" $JOB_NAME $status_file

check_exit_status "$java -Djava.io.tmpdir=$tempDir -Xmx15g -jar $gatk \
    -T PrintReads \
    -R $genome_fasta \
    -I $outDir/Aligned.sortedByCoord.out.dedupped.split.realigned.bam \
    -BQSR $outDir/recal_data.table \
    -o $outDir/Aligned.sortedByCoord.out.dedupped.split.realigned.recal.bam" $JOB_NAME $status_file

# VARIANT CALLING
check_exit_status "$java -Djava.io.tmpdir=$tempDir -Xmx15g -jar $gatk \
	-T HaplotypeCaller \
	-R $genome_fasta \
	-I $outDir/Aligned.sortedByCoord.out.dedupped.split.realigned.recal.bam \
	-dontUseSoftClippedBases \
	-stand_call_conf 20.0 \
	-o $outDir/${fastq_end1}_raw.vcf" $JOB_NAME $status_file

# TODO: VARIANT FILTERING
check_exit_status "$java -Djava.io.tmpdir=$tempDir -Xmx15g -jar $gatk \
	-T VariantFiltration \
	-R $genome_fasta \
	-V $outDir/${fastq_end1}_raw.vcf \
	-window 35 \
	-cluster 3 \
	-filterName FS \
	-filter "FS > 30.0" \
	-filterName QD \
	-filter "QD < 2.0" \
	-o $outDir/${fastq_end1}_filt.vcf" $JOB_NAME $status_file

# zip and index
# check_exit_status "$bgzip $outDir/${fastq_end1}_filt.vcf -c > $outDir/${fastq_end1}_filt.vcf.gz" $JOB_NAME $status_file
# check_exit_status "$tabix -p vcf $outDir/${fastq_end1}_filt.vcf.gz" $JOB_NAME $status_file


# Upload
aws s3 cp $outDir/ $output_address/ --recursive
##END_UPLOAD##


