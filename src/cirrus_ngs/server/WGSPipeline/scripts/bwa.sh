#!/bin/bash

input_dir=$1
sample_dir=$2
file1_name=$3
file2_name=$4
log_dir=$5

exec 1>>$log_dir/bwa.log
exec 2>>$log_dir/bwa.log

sample_dir="$sample_dir/bwa"

file1_name=${file1_name//.gz/}
file1_name=${file1_name//./.trim.}
file2_name=${file2_name//.gz/}
file2_name=${file2_name//./.trim.}

output=${file1_name//.*/.bam}

genome=/shared/workspace/software/genomes/Hsapiens/bwa/ucsc.hg19.fasta
bwa=/shared/workspace/software/bwa/bwa-0.7.12/bwa
samblaster=/shared/workspace/software/samblaster/samblaster
samtools=/shared/workspace/software/samtools/samtools-1.1/samtools

mkdir -p $sample_dir

if [ "$file2_name" == "NULL" ]; then
    if [ ! -f $sample_dir/$output ]; then
        echo "Performing single end alignment for $file1_name ..."

        $bwa mem -M -t 8 -R "@RG\tID:1\tPL:temp\tPU:tempID\tSM:tempPheno" -v 1 \
            $genome $input_dir/$file1_name | $samblaster | $samtools view \
            -Sb - > $sample_dir/$output

        echo "Finished single end alignment for $file1_name"
    else
        echo "Single end alignment has already been performed for $file1_name"
    fi
else
    if [ ! -f $sample_dir/$output ]; then
        echo "Performing paired end alignment for $file1_name and $file2_name ..."

        $bwa mem -M -t 8 -R "@RG\tID:1\tPL:temp\tPU:tempID\tSM:tempPheno" -v 1 \
            $genome $input_dir/$file1_name $input_dir/$file2_name | \
            $samblaster | $samtools view -Sb - > $sample_dir/$output

        echo "Finished paired end alignment for $file1_name and $file2_name"
    else
        echo "Paired end alignment has already been performed for $file1_name and $file2_name"
    fi
fi

##DEBUG##
echo "dir is $sample_dir"
echo `ls $sample_dir`
##ENDDEBUG## 

echo
