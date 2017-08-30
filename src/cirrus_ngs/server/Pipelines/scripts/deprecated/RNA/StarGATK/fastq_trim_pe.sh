#!/bin/sh

inputFile=$1

if [ ! -f $inputFile"_1.trim.fastq" ]; then
   java -jar /shared/workspace/software/Trimmomatic-0.30/trimmomatic-0.30.jar PE -threads 5 -phred33 -trimlog /shared/workspace/software/Trimmomatic-0.30/trimlog.log $inputFile"_1.fastq" $inputFile"_2.fastq" $inputFile"_1.trim.fastq" $inputFile"_1.unpaired.fastq" $inputFile"_2.trim.fastq" $inputFile"_2.unpaired.fastq" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

fi
