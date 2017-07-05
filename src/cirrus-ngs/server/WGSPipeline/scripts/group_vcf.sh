#!/bin/bash

outputFile=$1

# redirecting all output to a file
exec 1>>$2/$HOSTNAME"_groupVCF.o"
exec 2>>$2/$HOSTNAME"_groupVCF.o"

workspace=/scratch/workspace
rootdir=/shared/workspace/software
gatk=$rootdir/gatk/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar
bgzip=$rootdir/tabix-0.2.6/bgzip
tabix=$rootdir/tabix-0.2.6/tabix
genomeSeq=$rootdir/sequences/Hsapiens/ucsc.hg19.fasta
dbsnp138=$rootdir/variation/dbsnp_138.hg19.vcf

tmpDir=$workspace/tmp

# store arguments in a special array 
args=("$@")
# get number of elements 
ELEMENTS=${#args[@]}

for (( i=2;i<$ELEMENTS;i++)); do
    vcfFile="--variant $workspace/${args[${i}]} "
    fileList=$fileList$vcfFile
done

if [ ! -f $workspace/$outputFile.merged.vcf ]; then
   java -Xmx2g -Djava.io.tmpdir=$tmpDir \
	-jar $gatk \
	-T CombineGVCFs \
	-R $genomeSeq \
	$fileList \
	-o $workspace/$outputFile.merged.vcf
fi

if [ ! -f $workspace/$outputFile.g.vcf ]; then
   java -Xms454m -Xmx3181m -Djava.io.tmpdir=$tmpDir \
	-jar $gatk \
	-T GenotypeGVCFs \
	-R $genomeSeq \
	--variant $workspace/$outputFile.merged.vcf \
        -o $workspace/$outputFile.g.vcf \
        --dbsnp $dbsnp138
fi

if [ ! -f $workspace/$outputFile.g.vcf.gz ]; then
   $bgzip -c $workspace/$outputFile.g.vcf > $workspace/$outputFile.g.vcf.gz
fi

if [ ! -f $workspace/$outputFile.g.vcf.gz.tbi ]; then
   $tabix -p vcf $workspace/$outputFile.g.vcf.gz
fi


