#!/bin/bash

output_file=$2

perl /shared/workspace/RNASeqPipeline/kallisto_deseq_workflow/scripts/transcript2gene_kallisto_counts_Hsa.pl $output_file
