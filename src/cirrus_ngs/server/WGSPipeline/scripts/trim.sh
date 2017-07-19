#!/bin/bash

sample_dir=$1
upload_dir=$2
file1_name=$3
file2_name=$4
log_dir=$5


exec 1>>$log_dir/trim.log
exec 2>>$log_dir/trim.log


file1_name=${file1_name//.gz/}
file2_name=${file2_name//.gz/}

full_trim_log=$log_dir/trimlogs/

trim1_out=${file1_name//./.trim.}
unpair1_out=${file1_name//./.unpaired.}
trim2_out="${file2_name//./.trim.}"
unpair2_out=${file2_name//./.unpaired.}

mkdir -p $sample_dir
mkdir -p $full_trim_log

full_trim_log=$full_trim_log/${file1_name//.*/.trimlog}

if [ "$file2_name" != "NULL" ]; then
    if [[ ! -f $sample_dir/$trim1_out && ! -f $sample_dir/$trim2_out ]]; then
        echo "Performing paired end trimming on $file1_name and $file2_name ..."
        echo

        java -jar /shared/workspace/software/Trimmomatic-0.36/trimmomatic-0.36.jar \
            PE -threads 16 -phred33 -trimlog $full_trim_log $sample_dir/$file1_name \
            $sample_dir/$file2_name $sample_dir/$trim1_out $sample_dir/$unpair1_out \
            $sample_dir/$trim2_out $sample_dir/$unpair2_out LEADING:3 TRAILING:3 \
            SLIDINGWINDOW:4:15 MINLEN:36

        echo
        echo "Finished trimming $file1_name and $file2_name"
    else
        echo "Paired end trimming has already been done on $file1_name and $file2_name"
    fi
else
    if [ ! -f $sample_dir/$trim1_out ]; then
        echo "Performing single end trimming on $file1_name ..."
        echo

        java -jar /shared/workspace/software/Trimmomatic-0.36/trimmomatic-0.36.jar \
            SE -threads 16 -phred33 -trimlog /dev/stderr $sample_dir/$file1_name \
            $sample_dir/$trim1_out LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

        echo
        echo "Finished trimming $file1_name"
    else
        echo "Single end trimming has already been done on $file1_name"
    fi
fi

rm $sample_dir/$file1_name
rm $sample_dir/$file2_name

echo
echo "Uploading trimmed files to $upload_dir ..."
echo
sed "s/\r/\n/g" <<< `aws s3 sync $sample_dir $upload_dir`
echo
echo "Finished uploading trimmed files"

##DEBUG##
echo
echo "dir is $sample_dir"
echo `ls $sample_dir`
##ENDDEBUG##
echo
