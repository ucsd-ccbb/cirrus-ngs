import sys
import os
import subprocess
import PBSTracker
import YamlFileReader
import re

#ROOT_DIR = "/shared/workspace/WGSPipeline"
#PROJECT_DIR = "/scratch/{}"
ROOT_DIR = "/scratch"
SCRIPTS = "/shared/workspace/Pipelines/scripts/"
LOG_DIR = "/shared/workspace/logs/DNASeq/"
DATA_DIR = ""
CHROMOSOME_LIST = list(map(str, range(1,23))) + ["X", "Y", "M"]

##run the actual pipeline
def run_analysis(yaml_file):
    documents = YamlFileReader.parse_yaml_file(yaml_file)

    project_name = documents.get("project")
    analysis_steps = documents.get("analysis")
    output_address = documents.get("upload")
    sample_list = documents.get("sample")

    global LOG_DIR
    LOG_DIR += project_name

    if "fastqc" in analysis_steps:
        run_fastqc(project_name, sample_list, output_address)
    if not "notrim" in analysis_steps:
        run_trim(project_name, sample_list, output_address)
    if "bwa" in analysis_steps:
        run_bwa(project_name, sample_list, output_address)
    if "sort" in analysis_steps:
        run_sort(project_name, sample_list, output_address)
    if "dedup" in analysis_steps:
        run_dedup(project_name, sample_list, output_address)
    if "split" in analysis_steps:
        run_split(project_name, sample_list, output_address)
    if "postalignment" in analysis_steps:
        run_postalignment(project_name, sample_list, output_address)
    if "haplotype" in analysis_steps:
        run_gatkhaplotype(project_name, sample_list, output_address)
        



def run_fastqc(project_name, sample_list, output_address):
    subprocess_call_list = ["qsub", "-o", "/dev/null", "-e", "/dev/null", SCRIPTS + "fastqc.sh", project_name]
    for curr_sample_arguments in _sample_argument_generator(sample_list, output_address):
        argument_list = list(curr_sample_arguments)
        subprocess.call(subprocess_call_list + argument_list)

    PBSTracker.trackPBSQueue(1, "fastqc.sh")

    print("Finished running FastQC")

def run_trim(project_name, sample_list, output_address):
    subprocess_call_list = ["qsub", "-o", "/dev/null", "-e", "/dev/null", SCRIPTS + "trim.sh", project_name]
    for curr_sample_arguments in _sample_argument_generator(sample_list, output_address):
        argument_list = list(curr_sample_arguments)
        argument_list.append("4")
        argument_list.append("36")
        subprocess.call(subprocess_call_list + argument_list)

    PBSTracker.trackPBSQueue(1, "trim.sh")

    print("Finished running Trimmomatic")

def run_bwa(project_name, sample_list, output_address):
    subprocess_call_list = ["qsub", "-o", "/dev/null", "-e", "/dev/null", SCRIPTS + "bwa.sh", project_name]
    for curr_sample_arguments in _sample_argument_generator(sample_list, output_address):
        argument_list = list(curr_sample_arguments)
        argument_list[0] = ".trim" + argument_list[0]
        argument_list[4] = output_address
        argument_list[7] = "False"
        argument_list.append("8")
        subprocess.call(subprocess_call_list + argument_list)

    PBSTracker.trackPBSQueue(1, "bwa.sh")

    print("Finished running bwa")

def run_sort(project_name, sample_list, output_address):
    subprocess_call_list = ["qsub", "-o", "/dev/null", "-e", "/dev/null", SCRIPTS + "sort.sh", project_name]
    for curr_sample_arguments in _sample_argument_generator(sample_list, output_address):
        argument_list = list(curr_sample_arguments)
        argument_list[0] = ".bam"
        argument_list[4] = output_address
        argument_list[7] = "False"
        argument_list.append("4")
        subprocess.call(subprocess_call_list + argument_list)

    PBSTracker.trackPBSQueue(1, "sort.sh")

    print("Finished running sort")

def run_dedup(project_name, sample_list, output_address):
    subprocess_call_list = ["qsub", "-o", "/dev/null", "-e", "/dev/null", SCRIPTS + "dedup.sh", project_name]
    for curr_sample_arguments in _sample_argument_generator(sample_list, output_address):
        argument_list = list(curr_sample_arguments)
        argument_list[0] = ".sort.bam"
        argument_list[4] = output_address
        argument_list[7] = "False"
        argument_list.append("4")
        subprocess.call(subprocess_call_list + argument_list)

    PBSTracker.trackPBSQueue(1, "dedup.sh")

    print("Finished running dedup")

def run_split(project_name, sample_list, output_address):
    subprocess_call_list = ["qsub", "-o", "/dev/null", "-e", "/dev/null", SCRIPTS + "split.sh", project_name]
    for curr_sample_arguments in _sample_argument_generator(sample_list, output_address):
        argument_list = list(curr_sample_arguments)
        argument_list[0] = ".dedup.bam"
        argument_list[4] = output_address
        argument_list[7] = "False"
        argument_list.append("1")

        for chromosome in ["1", "2"]: #CHROMOSOME_LIST:
            subprocess.call(subprocess_call_list + argument_list + [chromosome])

    PBSTracker.trackPBSQueue(1, "split.sh")

    print("Finished running split")

def run_postalignment(project_name, sample_list, output_address):
    subprocess_call_list = ["qsub", "-o", "/dev/null", "-e", "/dev/null", SCRIPTS + "post.sh", project_name]
    for curr_sample_arguments in _sample_argument_generator(sample_list, output_address):
        argument_list = list(curr_sample_arguments)
        argument_list[0] = ".bam"
        argument_list[4] = output_address
        argument_list[7] = "False"
        argument_list.append("4")

        for chromosome in ["1", "2"]: #CHROMOSOME_LIST:
            subprocess.call(subprocess_call_list + argument_list + [chromosome])

    PBSTracker.trackPBSQueue(1, "post.sh")

    print("Finished running post")

def run_gatkhaplotype(project_name, sample_list, output_address):
    subprocess_call_list = ["qsub", "-o", "/dev/null", "-e", "/dev/null", SCRIPTS + "haplo.sh", project_name]
    for curr_sample_arguments in _sample_argument_generator(sample_list, output_address):
        argument_list = list(curr_sample_arguments)
        argument_list[4] = output_address
        argument_list[7] = "False"
        argument_list.append("4")

        for chromosome in ["1", "2"]: #CHROMOSOME_LIST:
            argument_list[0] = ".final.{}.bam".format(chromosome)
            subprocess.call(subprocess_call_list + argument_list + [chromosome])

    PBSTracker.trackPBSQueue(1, "haplo.sh")

    print("Finished running haplotype")

#returns tuple
#first element is file name without suffix
#second element is .fastq or .fq
#third element is boolean representing if file was zipped
def _separate_file_suffix(sample_file):
    #regex matches .fastq or .fq and then any extensions following them
    original_suffix = re.search("\.f(?:ast){0,1}q.*$", sample_file).group()
    file_prefix = sample_file.replace(original_suffix, "")
    is_zipped = ".gz" in original_suffix
    file_suffix = original_suffix.replace(".gz", "")

    return file_prefix, file_suffix, str(is_zipped)

#generator that yields tuple containing arguments for shell script
#yielded arguments are standard for every tool
#returned tuple:
#   file_suffix:    file extension (.fq or .fastq) without .gz
#   ROOT_DIR:       directory under which output will be saved
#   fastq_end1:     forward reads (or single end file)
#   fastq_end2:     reverse reads (or "NULL" if single end)
#   input_address:  from "download" value of each sample, is S3 address
#   output_address: S3 address where final analysis will be uploaded
#   LOG_DIR:        directory where logs will be stored
#   is_zipped:      str version of boolean, indicates if files to be downloaded are gzipped
def _sample_argument_generator(sample_list, output_address):
    for sample_pair in sample_list:
        curr_samples = [file_name.strip() for file_name in sample_pair.get("filename").split(",")]
        fastq_end1, file_suffix, is_zipped = _separate_file_suffix(curr_samples[0])
        
        #puts "NULL" in fastq_end2 if sample isn't paired end
        if len(curr_samples) > 1:
            fastq_end2 = _separate_file_suffix(curr_samples[1])[0]
        else:
            fastq_end2 = "NULL"

        yield (file_suffix, ROOT_DIR, fastq_end1, fastq_end2, 
                sample_pair.get("download"), output_address, LOG_DIR, is_zipped)


if __name__ == "__main__":
    run_analysis(sys.argv[1])
