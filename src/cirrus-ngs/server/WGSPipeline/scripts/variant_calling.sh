#!/bin/bash

normal_file=$1
tumor_file=$2
chromosome=$3

workspace=/scratch/workspace

# redirecting all output to a file
exec 1>>$workspace/"variant_call.o"
exec 2>>$workspace/"variant_call.e"

echo "$(date): running mutect..."
## MuTect
if [ ! -f $workspace/mutect.$chromosome.vcf.gz ]; then

   mkdir $workspace/temp/$chromosome
   mkdir $workspace/$normal_file"_vs_"$tumor_file

   java -Djava.io.tmpdir=$workspace/temp -Xmx4g -jar /shared/workspace/software/mutect/muTect-1.1.5.jar --analysis_type MuTect --reference_sequence /shared/shuling/resource/bwa/human_g1k_v37.fasta --cosmic /shared/shuling/resource/b37_cosmic_v54_120711.vcf --dbsnp /shared/shuling/resource/dbsnp_132_b37.leftAligned.vcf --input_file:normal $workspace/$normal_file/$normal_file.final.$chromosome.bam --input_file:tumor $workspace/$tumor_file/$tumor_file.final.$chromosome.bam -vcf $workspace/$normal_file"_vs_"$tumor_file/mutect.$chromosome.vcf.gz -L $chromosome
fi


