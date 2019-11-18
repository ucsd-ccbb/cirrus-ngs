#!/bin/bash

project_name=$1
workflow=$2
file_suffix=$3  #extension of input file, does not include .gz if present in input
root_dir=$4
fastq_end1=$5
fastq_end2=$6
input_address=$7    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$8   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$9
is_zipped=${10}    #either "True" or "False", indicates whether input is gzipped
num_threads=${11}       #number of threads

#logging
log_dir=$log_dir/$fastq_end1
mkdir -p $log_dir
log_file=$log_dir/'gatk_vc.log'
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$workflow/$fastq_end1
tempDir=$workspace/tmp

mkdir -p $tempDir

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME $status_file

##DOWNLOAD##
if [ ! -f $workspace/$fastq_end1$file_suffix ]
then
    aws s3 cp $input_address/$fastq_end1$file_suffix $workspace/ --quiet
    aws s3 cp $input_address/$fastq_end1$file_suffix.bai $workspace/ --quiet
fi
##END_DOWNLOAD##

# START THE PIPELINE

# MARK DUPLICATES
check_exit_status "$java -jar $picard_mark_duplicates \
	I=$workspace/$fastq_end1$file_suffix \
	O=$workspace/$fastq_end1.Aligned.sortedByCoord.out.dedupped.bam \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=SILENT \
	M=$workspace/$fastq_end1.output.metrics" $JOB_NAME $status_file

check_exit_status "$java -Djava.io.tmpdir=$tempDir -jar $picard_sort_sam \
    	I=$workspace/$fastq_end1.Aligned.sortedByCoord.out.dedupped.bam \
    	O=$workspace/$fastq_end1.Aligned.sorted.out.dedupped.bam \
        SORT_ORDER=coordinate \
	CREATE_INDEX=true" $JOB_NAME $status_file

check_exit_status "$java -Djava.io.tmpdir=$tempDir -Xmx15g -jar $gatk_v38 \
	-T SplitNCigarReads \
	-R $genome_fasta \
	-I $workspace/$fastq_end1.Aligned.sorted.out.dedupped.bam \
	-o $workspace/$fastq_end1.Aligned.sorted.out.dedupped.split.bam \
	-rf ReassignOneMappingQuality \
	-RMQF 255 \
	-RMQT 60 \
	-fixNDN \
	-U ALLOW_N_CIGAR_READS" $JOB_NAME $status_file


# INDEL REALIGNMENT
check_exit_status "$java -Djava.io.tmpdir=$tempDir -Xmx15g -jar $gatk_v38 \
	-T RealignerTargetCreator \
	-R $genome_fasta \
	-I $workspace/$fastq_end1.Aligned.sorted.out.dedupped.split.bam \
	-o $workspace/$fastq_end1.output.intervals \
	-known $mills \
	-known $G1000_snps" $JOB_NAME $status_file

check_exit_status "$java -Djava.io.tmpdir=$tempDir -Xmx15g -jar $gatk_v38 \
	-I $workspace/$fastq_end1.Aligned.sorted.out.dedupped.split.bam \
	-R $genome_fasta \
	-T IndelRealigner \
	-targetIntervals $workspace/$fastq_end1.output.intervals \
	-known $mills \
	-known $G1000_snps \
	-o $workspace/$fastq_end1.Aligned.sorted.out.dedupped.split.realigned.bam" $JOB_NAME $status_file

# BQSR
check_exit_status "$java -Djava.io.tmpdir=$tempDir -Xmx15g -jar $gatk_v38 \
    	-T BaseRecalibrator \
    	-R $genome_fasta \
    	-I $workspace/$fastq_end1.Aligned.sorted.out.dedupped.split.realigned.bam \
    	-knownSites $dbsnp \
    	-knownSites $mills \
    	-knownSites $G1000_snps \
    	-o $workspace/$fastq_end1.recal_data.table" $JOB_NAME $status_file

check_exit_status "$java -Djava.io.tmpdir=$tempDir -Xmx15g -jar $gatk_v38 \
    	-T PrintReads \
    	-R $genome_fasta \
    	-I $workspace/$fastq_end1.Aligned.sorted.out.dedupped.split.realigned.bam \
    	-BQSR $workspace/$fastq_end1.recal_data.table \
    	-o $workspace/$fastq_end1.Aligned.sorted.out.dedupped.split.realigned.bqsr.bam" $JOB_NAME $status_file

# VARIANT CALLING
check_exit_status "$java -Djava.io.tmpdir=$tempDir -Xmx15g -jar $gatk_v38 \
	-T HaplotypeCaller \
	-R $genome_fasta \
	-I $workspace/$fastq_end1.Aligned.sorted.out.dedupped.split.realigned.bqsr.bam \
	-dontUseSoftClippedBases \
	-U ALLOW_N_CIGAR_READS \
	-stand_call_conf 20.0 \
	--dbsnp $dbsnp \
	-out_mode EMIT_VARIANTS_ONLY \
	-rf BadCigar \
	-o $workspace/${fastq_end1}_raw.vcf \
	-nct $num_threads" $JOB_NAME $status_file

# VARIANT FILTRATION
$java -Djava.io.tmpdir=$tempDir -Xmx15g -jar $gatk_v38 \
	-T VariantFiltration \
	-R $genome_fasta \
	-V $workspace/${fastq_end1}_raw.vcf \
	-window 35 \
	-cluster 3 \
	-filterName FS -filter "FS > 30.0" \
	-filterName QD -filter "QD < 2.0" \
	-filterName Quality -filter "QUAL < 20" \
	--genotypeFilterExpression "DP < 5" --genotypeFilterName DP \
	--setFilteredGtToNocall \
	-o $workspace/${fastq_end1}_tmp.vcf

# SELECT VARIANTS
$java -Djava.io.tmpdir=$tempDir -Xmx15g -jar $gatk_v38 \
        -T SelectVariants \
        -R $genome_fasta \
        -V $workspace/${fastq_end1}_tmp.vcf \
        --excludeNonVariants \
	--excludeFiltered \
	-o $workspace/${fastq_end1}_filt.vcf

# zip and index
check_exit_status "$bgzip $workspace/${fastq_end1}_filt.vcf -c > $workspace/${fastq_end1}.vcf.gz" $JOB_NAME $status_file
check_exit_status "$tabix -p vcf $workspace/${fastq_end1}.vcf.gz" $JOB_NAME $status_file

##UPLOAD## 
aws s3 cp $workspace/ $output_address/ --recursive --quiet
##END_UPLOAD##


