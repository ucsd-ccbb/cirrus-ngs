#!/bin/bash

export R_LIBS="/shared/workspace/software/R-packages"

Rscript /shared/workspace/RNASeqPipeline/kallisto_deseq_workflow/scripts/RNA-seq_limma.R $1 $2
