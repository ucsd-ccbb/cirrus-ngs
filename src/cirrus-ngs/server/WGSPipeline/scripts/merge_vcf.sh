#!/bin/bash

export PERL5LIB=/shared/workspace/software/vcftools_0.1.12b/perl:$PERL5LIB
export PATH=/shared/workspace/software/vcftools_0.1.12b/perl:${PATH}
export PATH=/shared/workspace/software/vcftools_0.1.12b/bin:${PATH}
export PATH=/shared/workspace/software/tabix-0.2.6:${PATH}

workspace=/scratch/workspace
rootdir=/shared/workspace/software
bgzip=$rootdir/tabix-0.2.6/bgzip
tabix=$rootdir/tabix-0.2.6/tabix
vcfConcat=$rootdir/vcftools_0.1.12b/bin/vcf-concat
vcfSort=$rootdir/vcftools_0.1.12b/bin/vcf-sort

outputFile=$1
byRegion=$2

# redirecting all output to a file
exec 1>>$3/$HOSTNAME"_mergevcf.o"
exec 2>>$3/$HOSTNAME"_mergevcf.o"

tmpDir=$workspace/tmp

if [ ! -d $tmpDir ]; then
   mkdir $tmpDir
fi

# store arguments in a special array 
args=("$@")
# get number of elements 
ELEMENTS=${#args[@]}

for (( i=3;i<$ELEMENTS;i++)); do
    vcfFile="$workspace/${args[${i}]} "
    fileList=$fileList$vcfFile
done

if [ "$byRegion" == "true" ] && [ ! -f $workspace/$outputFile.raw.vcf.gz ]; then

   $vcfConcat $fileList | $bgzip -c > $workspace/$outputFile.vcf.gz

   $vcfSort -t $tmpDir $workspace/$outputFile.vcf.gz > $workspace/$outputFile.raw.vcf

   $bgzip -c $workspace/$outputFile.raw.vcf > $workspace/$outputFile.raw.vcf.gz
   
fi

if [ "$byRegion" == "false" ] && [ ! -f $workspace/$outputFile.final.vcf.gz ]; then
    
   $vcfConcat $fileList | $bgzip -c > $workspace/$outputFile.final.vcf.gz
   
   $tabix -p vcf $workspace/$outputFile.final.vcf.gz
fi
