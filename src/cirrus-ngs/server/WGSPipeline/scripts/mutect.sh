#!/bin/bash


rootdir=/scratch/workspace

mutect=$rootdir/share/java/mutect/muTect-1.1.5.jar
genomeSeq=$rootdir/genomes/Hsapiens/hg19/genome.fa
dbsnp138=$rootdir/genomes/Hsapiens/GRCh37/variation/dbsnp_138.vcf.gz
cosmicVCF=$rootdir/genomes/Hsapiens/GRCh37/variation/cosmic-v68-GRCh37.vcf.gz

tmpDir=/scratch/workspace/tmp

if [ ! -f $tumorBAM-mutect.vcf.gz ]; then
	java -Xms454m -Xmx4g -XX:+UseSerialGC -Djava.io.tmpdir=$tmpDir -jar $mutect -R $genomeSeq -T MuTect -U ALLOW_N_CIGAR_READS --downsample_to_coverage 10000 --read_filter NotPrimaryAlignment -I:tumor $tumorBAM.final.bam --tumor_sample_name syn3-tumor -I:normal $normalBAM.final.bam --normal_sample_name syn3-normal --fraction_contamination 0 --dbsnp $dbsnp138 --cosmic $cosmicVCF -L $normalBAM.final.bed --interval_set_rule INTERSECTION --vcf $tumorBAM-mutect.vcf.gz -o /dev/null
fi


