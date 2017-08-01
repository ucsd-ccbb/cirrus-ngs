#!/bin/sh

inputFile=$1

if [ ! -f $inputFile"_1.trim.fastq" ]; then
   java -jar /shared/workspace/software/Trimmomatic-0.30/trimmomatic-0.30.jar SE -threads 5 -phred33 -trimlog /shared/workspace/software/Trimmomatic-0.30/trimlog.log $inputFile".fastq" $inputFile".trim.fastq" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:27
fi

