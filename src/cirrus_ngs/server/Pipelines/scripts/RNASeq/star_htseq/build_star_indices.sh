#!/bin/bash

/shared/workspace/software/STAR/2.5.3a/bin/Linux_x86_64/STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /shared/workspace/software/references/Hsapiens/hg19/indices/STAR/ucsc/ --genomeFastaFiles /shared/workspace/software/references/Hsapiens/hg19/sequence/ucsc.hg19.fasta --sjdbGTFfile /shared/workspace/software/references/Hsapiens/hg19/annotations/gencode.v19.annotation.gtf --sjdbOverhang 74


#/shared/workspace/software/STAR/2.5.3a/bin/Linux_x86_64/STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /shared/workspace/software/references/Hsapiens/hg38/indices/STAR/ --genomeFastaFiles /shared/workspace/software/references/Hsapiens/hg38/sequence/GRCh38.p12.genome.fa --sjdbGTFfile /shared/workspace/software/references/Hsapiens/hg38/annotations/gencode.v29.annotation.gtf --sjdbOverhang 74
