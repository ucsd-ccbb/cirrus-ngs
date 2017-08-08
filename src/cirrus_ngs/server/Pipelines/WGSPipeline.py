import sys
import os
import subprocess
from util import PBSTracker
from util import YamlFileReader
import re
import yaml

ROOT_DIR = "/scratch"
SCRIPTS = "/shared/workspace/Pipelines/scripts/"
LOG_DIR = "/shared/workspace/logs/DNASeq/"
CHROMOSOME_LIST = list(map(str, range(1,23))) + ["X", "Y", "M"]

##run the actual pipeline
def run_analysis(yaml_file):
    documents = YamlFileReader.parse_yaml_file(yaml_file)

    project_name = documents.get("project")
    analysis_steps = documents.get("analysis")
    output_address = documents.get("upload")
    sample_list = documents.get("sample")

    group_list = {}

    for sample_pair in sample_list:
        curr_group = sample_pair.get("group")
        curr_sample = sample_pair.get("description")
        if group_list.get(curr_group, None):
            group_list.get(curr_group).append(curr_sample)
        else:
            group_list[curr_group] = [curr_sample]


    global LOG_DIR
    LOG_DIR += project_name

    config_file = open("/shared/workspace/Pipelines/tools.yaml", "r")
    tool_config_dict = yaml.load(config_file)

    wgs_specific_config_file = open("/shared/workspace/Pipelines/WGSTools.yaml", "r")
    wgs_config_dict = yaml.load(wgs_specific_config_file)

    for step in wgs_config_dict["steps"]:
        if step in analysis_steps:
            run_tool(tool_config_dict[step], wgs_config_dict[step], project_name, sample_list, output_address, group_list)



def run_tool(tool_config_dict, extra_bash_args, project_name, sample_list, output_address, group_list):
    subprocess_call_list = ["qsub", "-o", "/dev/null", "-e", "/dev/null", 
            SCRIPTS + tool_config_dict["script_name"] + ".sh", project_name]
    
    extra_bash_args = list(map(str, extra_bash_args))

    if tool_config_dict.get("by_group", False):
        for group_arguments in _by_group_argument_generator(group_list, output_address, tool_config_dict):
            subprocess.call(subprocess_call_list + group_arguments + extra_bash_args)

        PBSTracker.trackPBSQueue(1, tool_config_dict["script_name"])
        return

    for curr_sample_arguments in _sample_argument_generator(sample_list, output_address, tool_config_dict):
        if tool_config_dict["uses_chromosomes"]:
            original_suffix = curr_sample_arguments[0]
            for chromosome in CHROMOSOME_LIST: 
                curr_sample_arguments[0] = original_suffix.format(chromosome)
                subprocess.call(subprocess_call_list + curr_sample_arguments + 
                        extra_bash_args + [chromosome])
        else:
            subprocess.call(subprocess_call_list + curr_sample_arguments + extra_bash_args)

    PBSTracker.trackPBSQueue(1, tool_config_dict["script_name"])

#returns tuple
#first element is file name without suffix
#second element is .fastq or .fq
#third element is str boolean representing if file was zipped
def _separate_file_suffix(sample_file):
    #regex matches .fastq or .fq and then any extensions following them
    #user shouldn't have "fastq" or "fq" in their file names before the ext
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
def _sample_argument_generator(sample_list, output_address, config_dictionary):
    download_suffix = config_dictionary["download_suffix"]
    input_is_output = config_dictionary["input_is_output"]
    can_be_zipped = config_dictionary["can_be_zipped"]
    uses_chromosomes = config_dictionary["uses_chromosomes"]

    for sample_pair in sample_list:
        curr_samples = [file_name.strip() for file_name in sample_pair.get("filename").split(",")]
        fastq_end1, file_suffix, is_zipped = _separate_file_suffix(curr_samples[0])
        
        #puts "NULL" in fastq_end2 if sample isn't paired end
        if len(curr_samples) > 1:
            fastq_end2 = _separate_file_suffix(curr_samples[1])[0]
        else:
            fastq_end2 = "NULL"

        #used if different precursor than .fq or .fastq is needed
        if download_suffix:
            if uses_chromosomes:
                file_suffix = download_suffix
            else:
                file_suffix = download_suffix.format(file_suffix)

        #for tools later in pipeline the precursors are 
        #downloaded from the output address
        if input_is_output:
            input_address = output_address
        else:
            input_address = sample_pair.get("download")

        #some precursors are never zipped
        if not can_be_zipped:
            is_zipped = "False"


        yield [file_suffix, ROOT_DIR, fastq_end1, fastq_end2, 
                input_address, output_address, LOG_DIR, is_zipped]

def _by_group_argument_generator(group_list, output_address, config_dictionary):
    download_suffix = config_dictionary["download_suffix"]

    for group in group_list:
        samples = " ".join(group_list[group])
        yield [download_suffix, ROOT_DIR, group, "NULL", output_address,
                output_address, LOG_DIR, "False", samples]





if __name__ == "__main__":
    run_analysis(sys.argv[1])
