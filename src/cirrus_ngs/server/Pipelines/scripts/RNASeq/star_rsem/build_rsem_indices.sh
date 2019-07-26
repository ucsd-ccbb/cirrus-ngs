#!/bin/bash

export PERL5LIB="/shared/workspace/software/perl-Env/1.04/lib/"

## Human hg19
/shared/workspace/software/RSEM/1.3.0/rsem-prepare-reference --gtf /shared/workspace/software/references/Hsapiens/hg19/annotations/gencode.v19.annotation.gtf --star --star-path /shared/workspace/software/STAR/2.5.3a/bin/Linux_x86_64 -p 8 /shared/workspace/software/references/Hsapiens/hg19/sequence/GRCh37.p13.genome.fa /shared/workspace/software/references/Hsapiens/hg19/indices/RSEM/GRCh37/human

## Human hg38
#/shared/workspace/software/RSEM/1.3.0/rsem-prepare-reference --gtf /shared/workspace/software/references/Hsapiens/hg38/annotations/gencode.v29.annotation.gtf --star --star-path /shared/workspace/software/STAR/2.5.3a/bin/Linux_x86_64 -p 8 /shared/workspace/software/references/Hsapiens/hg38/sequence/GRCh38.p12.genome.fa /shared/workspace/software/references/Hsapiens/hg38/indices/STAR/human

#/shared/workspace/software/RSEM/1.3.0/rsem-prepare-reference --gtf /shared/workspace/software/references/Hsapiens/hg38/annotations/gencode.v29.annotation.gtf --star --star-path /shared/workspace/software/STAR/2.5.3a/bin/Linux_x86_64 -p 8 /shared/workspace/software/references/Hsapiens/hg38/sequence/GRCh38.p12.genome.fa /shared/workspace/software/references/Hsapiens/hg38/indices/RSEM

## Mouse M19
#/shared/workspace/software/RSEM/1.3.0/rsem-prepare-reference --gtf /shared/workspace/software/references/Mmusculus/M19/annotations/gencode.vM19.annotation.gtf --star --star-path /shared/workspace/software/STAR/2.5.3a/bin/Linux_x86_64 -p 8 /shared/workspace/software/references/Mmusculus/M19/sequence/GRCm38.p6.genome.fa /shared/workspace/software/references/Mmusculus/M19/indices/STAR/mouse

## Mouse GRCm38
#/shared/workspace/software/RSEM/1.3.0/rsem-prepare-reference --gtf /shared/workspace/software/references/Mmusculus/GRCm38/annotations/Mus_musculus.GRCm38.68.gtf --star --star-path /shared/workspace/software/STAR/2.5.3a/bin/Linux_x86_64 -p 8 /shared/workspace/software/references/Mmusculus/GRCm38/sequence/GRCm38_68.fa /shared/workspace/software/references/Mmusculus/GRCm38/indices/STAR/mouse

## Rat v6
#/shared/workspace/software/RSEM/1.3.0/rsem-prepare-reference --gtf /shared/workspace/software/references/Rat/Rnor6/annotations/Rnor_6.0.94.chr.gtf --star --star-path /shared/workspace/software/STAR/2.5.3a/bin/Linux_x86_64 -p 8 /shared/workspace/software/references/Rat/Rnor6/sequence/Rnor_6.0_94.fa /shared/workspace/software/references/Rat/Rnor6/indices/STAR/rat

## caenorhabditis
#/shared/workspace/software/RSEM/1.3.0/rsem-prepare-reference --gtf /shared/workspace/software/references/caenorhabditis/WBcel235/annotations/Caenorhabditis_elegans.WBcel235.96.gtf --star --star-path /shared/workspace/software/STAR/2.5.3a/bin/Linux_x86_64 -p 8 /shared/workspace/software/references/caenorhabditis/WBcel235/sequence/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa /shared/workspace/software/references/caenorhabditis/WBcel235/indices/STAR/caenorhabditis
