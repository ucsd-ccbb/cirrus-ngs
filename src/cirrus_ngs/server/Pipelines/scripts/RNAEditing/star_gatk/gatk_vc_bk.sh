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
log_file=$log_dir/'gatk_filt.log'
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
    aws s3 cp $input_address/$fastq_end1.Aligned.sorted.out.dedupped.split.realigned.bqsr.bam $workspace/ --quiet
    aws s3 cp $input_address/$fastq_end1.Aligned.sorted.out.dedupped.split.realigned.bqsr.bai $workspace/ --quiet
fi
##END_DOWNLOAD##

### filter $workspace/${fastq_end1}_filt.vcf for FILTER=="PASS"
check_exit_status  "awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' $workspace/${fastq_end1}_filt.vcf > $workspace/${fastq_end1}_filt.vcf.pass" $JOB_NAME $status_file

# echo "Performing SNPir filtering steps, except for RNA editing site removal";
# convert vcf format into custom SNPiR format and filter variants with quality <20
check_exit_status "$SNPiR/convertVCF.sh $workspace/${fastq_end1}_filt.vcf.pass $workspace/${fastq_end1}_raw_variants.txt 20" $JOB_NAME $status_file
# filt for total depth of at least 5 reads and 2 reads with alternate allele
check_exit_status "awk -F '\t|,' '{if (($3 >= 5) && ($4 >= 2)) {print}}' $workspace/${fastq_end1}_raw_variants.txt > $workspace/${fastq_end1}_raw_variants.mincov.txt" $JOB_NAME $status_file

# filter mismatches at read ends
# note: add the -illumina option if your reads are in Illumina 1.3+ quality format
check_exit_status "$SNPiR/filter_mismatch_first6bp.pl \
	-infile $workspace/${fastq_end1}_raw_variants.mincov.txt \
	-outfile $workspace/${fastq_end1}_raw_variants.rmhex.txt \
	-bamfile $workspace/${fastq_end1}.Aligned.sorted.out.dedupped.split.realigned.bqsr.bam" $JOB_NAME $status_file

# # alu sites
check_exit_status "awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' $workspace/${fastq_end1}_raw_variants.rmhex.txt | \
	$intersectBed \
	-wa \
	-header \
	-a stdin \
	-b $hg19_AluRegions \
	> $workspace/${fastq_end1}_alu_sites.txt" $JOB_NAME $status_file 

# # alu sites from vars
check_exit_status "$intersectBed \
	-wa \
	-header \
	-a $workspace/${fastq_end1}_filt.vcf.pass \
	-b $workspace/${fastq_end1}_alu_sites.txt \
	> $workspace/${fastq_end1}_alu_rnaseq_vars.vcf" $JOB_NAME $status_file

# # alu RNA edit sites 
echo "Filtering for Alu RNA editing sites with RADAR and DARNED"
check_exit_status "$intersectBed \
	-wa \
	-a $workspace/${fastq_end1}_alu_sites.txt \
	-b $hg19_rnaEditDB \
	| uniq \
	> $workspace/${fastq_end1}_alu_rnaedit_sites.txt" $JOB_NAME $status_file

# # retrieve alu RNA edit sites from vars
check_exit_status "awk '{OFS="\t";print $1,$2,$2}' $workspace/${fastq_end1}_alu_rnaedit_sites.txt | \
	$intersectBed \
	-wa \
	-header \
	-a $workspace/${fastq_end1}_filt.vcf.pass \
	-b stdin \
	| uniq \
	> $workspace/${fastq_end1}_alu_rnaedit_sites.vcf" $JOB_NAME $status_file


# Continue on to nonalu sites
# filter out alu regions
check_exit_status "awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' $workspace/${fastq_end1}_raw_variants.rmhex.txt | \
	$intersectBed \
	-v \
	-a stdin \
	-b $hg19_AluRegions \
	> $workspace/${fastq_end1}_raw_variants.rmhex.nonalu.txt" $JOB_NAME $status_file

# filter variants in simple repeat regions
check_exit_status "awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' $workspace/${fastq_end1}_raw_variants.rmhex.nonalu.txt | \
	$intersectBed \
	-v \
	-a stdin \
	-b $hg19_RepeatMasker | \
	cut -f1,3-7 > $workspace/${fastq_end1}_raw_variants.rmhex.nonalu.rmsk.txt" $JOB_NAME $status_file

# filter intronic sites that are within 4bp of splicing junctions
# make sure your gene annotation file is in UCSC text format and sorted by chromosome and 
# transcript start position
check_exit_status "$SNPiR/filter_intron_near_splicejuncts.pl \
	-infile $workspace/${fastq_end1}_raw_variants.rmhex.nonalu.rmsk.txt \
	-outfile $workspace/${fastq_end1}_raw_variants.rmhex.nonalu.rmsk.rmintron.txt \
	-genefile $hg19_snpir_gene_annotation" $JOB_NAME $status_file

# filter variants in homopolymers
check_exit_status "$SNPiR/filter_homopolymer_nucleotides.pl \
	-infile $workspace/${fastq_end1}_raw_variants.rmhex.nonalu.rmsk.rmintron.txt \
	-outfile $workspace/${fastq_end1}_raw_variants.rmhex.nonalu.rmsk.rmintron.rmhom.txt \
	-refgenome $hg19_fasta" $JOB_NAME $status_file

# filter variants that were caused by mismapped reads
# this may take a while depending on the number of variants to screen and the size of the reference genome
# note: add the -illumina option if your reads are in Illumina 1.3+ quality format
check_exit_status "$SNPiR/BLAT_candidates.pl \
	-infile $workspace/${fastq_end1}_raw_variants.rmhex.nonalu.rmsk.rmintron.rmhom.txt \
	-outfile $workspace/${fastq_end1}_raw_variants.rmhex.nonalu.rmsk.rmintron.rmhom.rmblat.txt \
	-bamfile $workspace/${fastq_end1}.Aligned.sorted.out.dedupped.split.realigned.bqsr.bam \
	-refgenome $hg19_fasta" $JOB_NAME $status_file

##############################################################################################################

# filter for RNA Editing sites 

check_exit_status "awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' $workspace/${fastq_end1}_raw_variants.rmhex.nonalu.rmsk.rmintron.rmhom.rmblat.txt | \
	$intersectBed \
	-wa \
	-a stdin \
	-b $hg19_rnaEditDB \
	> $workspace/${fastq_end1}_nonalu_rnaedit_sites.txt" $JOB_NAME $status_file 

check_exit_status "awk '{OFS="\t";print $1,$2,$2}' $workspace/${fastq_end1}_raw_variants.rmhex.nonalu.rmsk.rmintron.rmhom.rmblat.txt | \
	$intersectBed \
	-wa \
	-header \
	-a $workspace/${fastq_end1}_filt.vcf.pass \
	-b stdin \
	| uniq \
	> $workspace/${fastq_end1}_nonalu_rnaseq_vars.vcf" $JOB_NAME $status_file

# Use RNA Editing site to intersect VCF file, annotate
# echo "Selecting RNA Editing site from VCF"
check_exit_status "$intersectBed \
	-wa \
	-header \
	-a $workspace/${fastq_end1}_filt.vcf.pass \
	-b $workspace/${fastq_end1}_nonalu_rnaedit_sites.txt \
	| uniq \
	> $workspace/${fastq_end1}_nonalu_rnaedit_sites.vcf" $JOB_NAME $status_file

##UPLOAD## 
aws s3 cp $workspace/ $output_address/ --recursive --quiet
##END_UPLOAD##

rm -r $workspace/

