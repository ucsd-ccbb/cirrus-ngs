1. Convert .sra to .csfasta/.qual with abi-dump
   Example: java -jar fileselector.jar /data/workspace/miRNASeq/shellScript/csconvert.sh /data/data_archive/smallRnaSeq/SRP012016_Han_Gastric_cancer_miRNAseq/SRA_folder/smrnaSeq/ .sra

2. Remove space from readID in .csfasta/.qual files with CSFastaConverter.jar
   Exmaple: java -jar CSFastaConverter.jar /data/data_archive/smallRnaSeq/SRP012016_Han_Gastric_cancer_miRNAseq/SRA_folder/smrnaSeq/ .csfasta
   Exmaple: java -jar CSFastaConverter.jar /data/data_archive/smallRnaSeq/SRP012016_Han_Gastric_cancer_miRNAseq/SRA_folder/smrnaSeq/ .qual

3. Convert .csfasta/.qual to .xsq with xsqconvert.sh
   Example: java -jar fileselector.jar /data/workspace/miRNASeq/shellScript/xsqconvert.sh /data/data_archive/smallRnaSeq/SRP012016_Han_Gastric_cancer_miRNAseq/SRA_folder/smrnaSeq/folder/ .csfasta

4. Run alignment with Lifescope



