__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import os
import re
import sys
import datetime

## parse software.config file and save the reference genomes and tools into dictionary.
def parse(workspace):
    reference_list = {}
    tool_list = {}

    config_file = workspace.replace("notebooks/cirrus-ngs", "src/cirrus_ngs") + "/server/Pipelines/config/software.conf"

    with open(config_file, 'r+') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("export"):
                config_key = line[len("export") + 1: line.find("=")]
                config_value = line[line.find("=") + 1: -1]

                if config_value.find("reference_dir") > -1:
                    reference_value = config_value[config_value.rfind("/") + 1: -1]
                    reference_list.update({config_key:reference_value})

                if config_value.find("tool_dir") > -1:
                    fields = re.split(r'/+', config_value)
                    if len(fields) > 3:
                        tool_value = fields[2]
                        tool_list.update({config_key:tool_value})

    return reference_list, tool_list


## according to pipeline, workflow and genome to print the information.
def print_software_info(pipeline, workflow, genome, reference_list, tool_list):
    print("#Primary analysis details")
    print("#Author: Guorong Xu")

    now = datetime.datetime.now()
    print("#Date: " + now.strftime("%Y-%m-%d %H:%M:%S") + "\n")

    print("#Reference used:")
    if (pipeline == "SmallRNASeq"):
        print("Reference genome: hairpin_" + genome)
    else:
        print("Reference genome: " + reference_list.get(genome + "_fasta"))
        print("Annotation: " + reference_list.get(genome + "_gtf") + "\n")

    print("#Tools used:")
    print("FASTQC: " + tool_list.get("fastqc"))
    print("Trimmomatic: " + tool_list.get("trimmomatic"))
    print("samtools: " + tool_list.get("samtools"))
    print("MultiQC: v1.3")

    ## RNASeq pipeline
    if (pipeline == "RNASeq"):
        if (workflow == "star_rsem"):
            print("STAR: " + tool_list.get("STAR"))
            print("RSEM: " + tool_list.get("rsem"))
        if (workflow == "star_gatk"):
            print("STAR: " + tool_list.get("STAR"))
            print("GATK: " + tool_list.get("gatk"))
        if (workflow == "star_htseq"):
            print("STAR: " + tool_list.get("STAR"))
            print("htseq: v0.9.1")
        if (workflow == "kallisto"):
            print("kallisto: " + tool_list.get("kallisto"))

    ## DNASeq pipeline
    if (pipeline == "DNASeq"):
        print("bwa: " + tool_list.get("bwa"))
        print("GATK: " + tool_list.get("gatk"))
        print("sambamba: " + tool_list.get("sambamba"))
        print("picard: " + tool_list.get("picard_add_or_replace_read_groups"))
        print("bedtools: " + tool_list.get("bedtools"))
        if (workflow == "bwa_mutect"):
            print("Mutect: mutect2")
        if (workflow == "bwa_sv"):
            print("lumpy: v0.2.14a-2")
            print("manta: " + tool_list.get("manta"))
            print("sv2: v1.4.3 ")

    ## RNAEditing pipeline
    if (pipeline == "RNAEditing"):
        print("bwa: " + tool_list.get("bwa"))
        print("GATK: " + tool_list.get("gatk"))
        print("sambamba: " + tool_list.get("sambamba"))
        print("picard: " + tool_list.get("picard_add_or_replace_read_groups"))
        print("bcftools: " + tool_list.get("bcftools"))
        print("tabix: " + tool_list.get("tabix"))
        print("bgzip: " + tool_list.get("bgzip"))


    ## ChipSeq pipeline
    if (pipeline == "ChiPSeq"):
        if (workflow == "homer"):
            print("bowtie: " + tool_list.get("bowtie"))
            print("homer: " + tool_list.get("make_tag_directory"))

    ## miRNASeq pipeline
    if (pipeline == "SmallRNASeq"):
        if (workflow == "bowtie2"):
            print("bowtie2: " + tool_list.get("bowtie2"))
            print("cutadapt: v1.14")

    print("\n\n")
