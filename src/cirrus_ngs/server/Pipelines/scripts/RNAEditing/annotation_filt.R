# Filter and merge Alu and Nonalu RNA edit oncotator files

library(data.table)
library(plyr, lib="/shared/workspace/software/anaconda2/lib/R/library")

args <- commandArgs(TRUE)
out.maf <- args[1]
known.alu <- args[2]
known.nonalu <- args[3]

filterRNAEdit <- function(rnaedit.table, novelty){
  rnaedit.table$Tumor_Sample_Barcode <- rnaedit.table$Matched_Norm_Sample_Barcode
  maf <- subset(rnaedit.table, (is.na(`1000gp3_AF`) | `1000gp3_AF` == ""))
  maf <- subset(maf, (is.na(ExAC_AF) | ExAC_AF == ""))
  maf <- subset(maf, !grepl("pseudogene", maf$gene_type))
  maf <- subset(maf, dbSNP_RS == "")
  maf$Novelty <- novelty
  maf <- merge(maf, strand.table, by = "Hugo_Symbol")
  maf$mutsig <- paste0(maf$Reference_Allele, ">", maf$Tumor_Seq_Allele2, maf$strand)
  maf <- subset(maf, mutsig %in%  c("A>G+", "T>C-"))
  maf$t_alt_count <- as.numeric(maf$t_alt_count)
  maf$t_ref_count <- as.numeric(maf$t_ref_count)
  maf$VAF <- maf$t_alt_count/(maf$t_alt_count + maf$t_ref_count)
  maf <- subset(maf, (t_alt_count >= 2 & (t_alt_count + t_ref_count >= 10)) & VAF > 0.1 & VAF < 0.95)
  maf <- unique(maf)
  maf
}
cols <- c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", "Variant_Classification", "Variant_Type", "Reference_Allele",
                "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "dbSNP_RS", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Protein_Change",
                "t_alt_count", "t_ref_count", "dbNSFP_CADD_phred", "ExAC_AF", "1000gp3_AF", "gene_type")
# Strand Data
strand.table <- unique(fread("/shared/workspace/software/references/Hsapiens/hg19/annotations/gene_strand_hg19.txt", stringsAsFactors = FALSE, sep = "\t", header=FALSE))
names(strand.table) <- c("Hugo_Symbol", "strand")

# RNA Edit data
# Known Alu
known.alu <- fread(known.alu, stringsAsFactors=FALSE, sep="\t", select=cols, skip="Hugo_Symbol")
known.alu  <- filterRNAEdit(known.alu, "Known Alu")
# Known Nonalu
known.nonalu <- fread(known.nonalu, stringsAsFactors=FALSE, sep="\t", select=cols, skip="Hugo_Symbol")
known.nonalu <- filterRNAEdit(known.nonalu, "Known Nonalu")

# Novel
if(!is.na(args[4])){
	# Novel Alu
	novel.alu <- args[4]
	novel.alu <- fread(novel.alu, stringsAsFactors=FALSE, sep="\t", select=cols, skip="Hugo_Symbol")
	novel.alu  <- filterRNAEdit(novel.alu, "Known Alu")
	# Novel Nonalu
	novel.nonalu <- args[5]
	novel.nonalu <- fread(novel.nonalu, stringsAsFactors=FALSE, sep="\t", select=cols, skip="Hugo_Symbol")
	novel.nonalu <- filterRNAEdit(novel.nonalu, "Known Nonalu")
	maf <- do.call(rbind, list(known.alu, known.nonalu, novel.alu, novel.nonalu))
} else {maf <- do.call(rbind, list(known.alu, known.nonalu))}

maf$Tumor_Sample_Barcode <- maf$Matched_Norm_Sample_Barcode
write.table(maf, file=out.maf, quote = FALSE, sep = "\t", row.names = FALSE)
