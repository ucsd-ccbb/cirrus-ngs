#!/bin/bash

project_name=$1
file_suffix=$2
file_path=$3
file1_name=$4
file2_name=$5

#redirect stdout and stderr to a log file
exec 1>>~/$project_name/bwa.log
exec 2>>~/$project_name/bwa.log

output=${file1_name//.*/.bam}

genome=/shared/workspace/software/genomes/Hsapiens/bwa/ucsc.hg19.fasta
bwa=/shared/workspace/software/bwa/bwa-0.7.12/bwa
samblaster=/shared/workspace/software/samblaster/samblaster
samtools=/shared/workspace/software/samtools/samtools-1.1/samtools


if [ "$file2_name" == "NULL" ]; then
    if [ ! -f $sample_dir/$output ]; then

        #if file2_name isn't NULL performs single end alignment

        echo "Performing single end alignment for $file1_name ..."

        $bwa mem -M -t 8 -R "@RG\tID:1\tPL:temp\tPU:tempID\tSM:tempPheno" -v 1 \
            $genome $sample_dir/$file1_name | $samblaster | $samtools view \
            -Sb - > $sample_dir/$output

        echo
        echo "Finished single end alignment for $file1_name"
    else
        echo "Single end alignment has already been performed for $file1_name"
    fi
else
    if [ ! -f $sample_dir/$output ]; then

        #if file2_name not NULL it does paired end alignment
        
        echo "Performing paired end alignment for $file1_name and $file2_name ..."

        $bwa mem -M -t 8 -R "@RG\tID:1\tPL:temp\tPU:tempID\tSM:tempPheno" -v 1 \
            $genome $sample_dir/$file1_name $sample_dir/$file2_name | \
            $samblaster | $samtools view -Sb - > $sample_dir/$output

        echo
        echo "Finished paired end alignment for $file1_name and $file2_name"
    else
        echo "Paired end alignment has already been performed for $file1_name and $file2_name"
    fi
fi

#the .fastq and .fq files will no longer be needed
rm $sample_dir/*.f*q

#uploads alignment results to s3 

echo
echo "Uploading aligned file to $upload_dir"
echo

#this is done in order to remove the carriage returns in output
sed "s/\r/\n/g" <<< `aws s3 sync --exclude *.f*q $sample_dir $upload_dir`

echo
echo "Finished uploading aligned file"


##DEBUG##
echo
echo "dir is $sample_dir"
echo `ls $sample_dir`
##ENDDEBUG## 

echo
