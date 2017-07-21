#!/bin/bash

project_name=$1
file_suffix=$2
file_path=$3
file1_name=$4
file2_name=$5

software_dir=/shared/workspace/software
workspace=$file_path/$project_name
trimmo_path=$software_dir/Trimmomatic-0.36/trimmomatic-0.36.jar

echo "This is in trim.sh "
echo $file_suffix
echo $file1_name
echo $file2_name

if [ ! -d $workspace/$file1_name ]; then
   mkdir -p $workspace/$file1_name
fi

# redirecting all output to a file
exec 1>>$workspace/'trim.log'
exec 2>>$workspace/'trim.log'

if [ "$file2_name" == "NULL" ];
then
    echo “Running trimmomatic on $file1_name...”
    java -jar $trimmo_path SE -threads 5 -phred33 -trimlog $workspace/'full_trimlogs' \
    $workspace/$file1_name/$file1_name$file_suffix \
    $workspace/$file1_name/$file1_name".trim"$file_suffix \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:27

    echo "Finished trimming $file1_name"

else
    echo “Running trimmomatic on $file1_name and $file2_name...”
    java -jar $trimmo_path SE -threads 5 -phred33 -trimlog $workspace/'full_trimlogs' \
    $workspace/$file1_name/$file1_name$file_suffix \
    $workspace/$file1_name/$file1_name".trim"$file_suffix \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:27

    java -jar $trimmo_path SE -threads 5 -phred33 -trimlog $workspace/'full_trimlogs' \
    $workspace/$file2_name/$file2_name$file_suffix \
    $workspace/$file2_name/$file2_name".trim"$file_suffix \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:27

    echo "Finished trimming $file1_name and $file2_name"
fi
