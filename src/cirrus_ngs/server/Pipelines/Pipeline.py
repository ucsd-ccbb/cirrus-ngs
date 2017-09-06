import sys
import os
import subprocess
from util import PBSTracker
#from util import YamlFileReader
import re
import yaml

ROOT_DIR = "/scratch"
SCRIPTS = "/shared/workspace/Pipelines/scripts/"
CHROMOSOME_LIST = list(map(str, range(1,23))) + ["X", "Y", "M"]

##run the actual pipeline
def run_analysis(yaml_file, log_dir, pipeline_config_file):
    #reads in the project's yaml file, YamlFileReader deprecated
    yaml_file_stream = open(yaml_file)
    documents = yaml.load(yaml_file_stream)
    yaml_file_stream.close()

    #from the sample agnostic region of the yaml file
    pipeline = documents.get("pipeline")
    project_name = documents.get("project")
    workflow = documents.get("workflow")
    analysis_steps = documents.get("analysis")
    output_address = documents.get("upload")
    input_address = documents.get("download")
    sample_list = documents.get("sample")

    os.environ["style"] = documents.get("style")
    os.environ["genome"] = documents.get("genome")

    if os.environ["style"] == "histone":
        os.environ["style_ext"] = ".regions.txt"
    elif os.environ["style"] == "factor":
        os.environ["style_ext"] = ".peaks.txt"


    #used for pair-based analysis
    #dictionary with key=normal_sample, value=tumor_sample
    pair_list = documents.get("pairs")

    os.environ["pairs_exist"] = str(not len(pair_list) == 0)


    #used for group-based analysis
    #dictionary with key=group, value=list of tuples from _separate_file_suffix
    group_list = {}

    for sample_pair in sample_list:
        curr_group = sample_pair.get("group")
        curr_sample = sample_pair.get("filename").split(",")[0].strip()
        curr_sample_tuple = _separate_file_suffix(curr_sample)
        if group_list.get(curr_group, None):
            group_list.get(curr_group).append(curr_sample_tuple)
        else:
            group_list[curr_group] = [curr_sample_tuple]


    #tools.yaml contains general configuration for all shell scripts
    #see comments in tools.yaml for more information
    global_config_file = open("/shared/workspace/Pipelines/config/tools.yaml", "r")
    global_config_dict = yaml.load(global_config_file)

    #configuration file for current pipeline
    #contains order of steps and extra arguments to each shell script
    specific_config_file = open("/shared/workspace/Pipelines/config/{}/{}".format(pipeline, pipeline_config_file), "r")
    specific_config_dict = yaml.load(specific_config_file)

    #steps list enforces order of possible steps in pipeline
    for step in specific_config_dict["steps"]:
        if step in analysis_steps:
            run_tool(global_config_dict[step], specific_config_dict[step], 
                    project_name, workflow, sample_list, input_address, output_address, group_list, pair_list, log_dir)

#function to run any tool
#arguments:
#   tool_config_dict: dictionary from parsing tools.yaml file
#   extra_bash_args: list from current step's value in current pipeline's yaml file
#   project_name: the project_name from the project's yaml file
#   sample_list: a list of dictionaries from parsing the sample field of the project's yaml file
#   input_address: the root s3 bucket that contains the original fastq/fq files for project
#   output_address: the root s3 bucket for uploading tool output
#   group_list: a dictionary as described in function run_analysis' group_list variable
#   pair_list: a dictionary as described in function run_analysis' pair_list variable
#   log_dir: root directory where logs will be stored
def run_tool(tool_config_dict, extra_bash_args, project_name, workflow, sample_list, input_address, output_address, group_list, pair_list, log_dir):
    #num_threads used to request nodes with a specific amount of threads for step
    if len(extra_bash_args) > 0:
        num_threads = extra_bash_args[0]
    else:
        num_threads = 1

    #contains all qsub flags and the name of the shell script for this step
    subprocess_call_list = ["qsub", "-V", "-o", "/dev/null", "-e", "/dev/null", "-pe", "smp", str(num_threads),
            SCRIPTS + tool_config_dict["script_path"] + ".sh"]
    
    #extra arguments to the current shell script outside of the 9 necessary ones (see shell script format documentation)
    #if not empty, first extra argument must be a number representing number of threads used by script
    extra_bash_args = list(map(str, extra_bash_args))

    #runs tool on every sample in project
    if tool_config_dict.get("all_samples", False):
        all_sample_arguments = _by_all_sample_argument_generator(project_name, workflow, sample_list, input_address, output_address, tool_config_dict, log_dir)
        if tool_config_dict["uses_chromosomes"]:
            original_suffix = all_sample_arguments[1]
            for chromosome in CHROMOSOME_LIST:
                all_sample_arguments[1] = original_suffix.format(chromosome)
                subprocess.call(subprocess_call_list + all_sample_arguments + extra_bash_args)
        else:
            subprocess.call(subprocess_call_list + all_sample_arguments + extra_bash_args)
        return

    if tool_config_dict.get("by_pair", None) and os.environ["pairs_exist"] == "True":
        for pair_arguments in _by_pair_argument_generator(project_name, workflow, group_list, pair_list, input_address, output_address, tool_config_dict, log_dir):
            if tool_config_dict["uses_chromosomes"]:
                original_suffix = pair_arguments[1]
                for chromosome in CHROMOSOME_LIST:
                    pair_arguments[1] = original_suffix.format(chromosome)
                    subprocess.call(subprocess_call_list + pair_arguments + extra_bash_args + [chromosome])
            else:
                subprocess.call(subprocess_call_list + pair_arguments + extra_bash_args)
        
        PBSTracker.trackPBSQueue(1, tool_config_dict["script_path"])
        return

    #runsn tool on samples in each group
    if tool_config_dict.get("by_group", False):
        for group_arguments in _by_group_argument_generator(project_name, workflow, group_list, input_address, output_address, tool_config_dict, log_dir):
            if tool_config_dict["uses_chromosomes"]:
                original_suffix = group_arguments[1]
                for chromosome in CHROMOSOME_LIST:
                    group_arguments[1] = original_suffix.format(chromosome)
                    subprocess.call(subprocess_call_list + group_arguments + 
                            extra_bash_args + [chromosome])
            else:
                subprocess.call(subprocess_call_list + group_arguments + extra_bash_args)

        PBSTracker.trackPBSQueue(1, tool_config_dict["script_path"])
        return

    for curr_sample_arguments in _sample_argument_generator(project_name, workflow, sample_list, input_address, output_address, tool_config_dict, log_dir):
        #for tools that run on each chromosome the file suffix has the current chrom number added to it
        if tool_config_dict["uses_chromosomes"]:
            original_suffix = curr_sample_arguments[1]
            for chromosome in CHROMOSOME_LIST:
                curr_sample_arguments[1] = original_suffix.format(chromosome)
                subprocess.call(subprocess_call_list + curr_sample_arguments + 
                        extra_bash_args + [chromosome])
        else:
            subprocess.call(subprocess_call_list + curr_sample_arguments + extra_bash_args)

    PBSTracker.trackPBSQueue(1, tool_config_dict["script_path"])

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
def _sample_argument_generator(project_name, workflow, sample_list, input_address, output_address, config_dictionary, log_dir):
    download_suffix = config_dictionary["download_suffix"]
    input_is_output = config_dictionary["input_is_output"]
    can_be_zipped = config_dictionary["can_be_zipped"]
    uses_chromosomes = config_dictionary["uses_chromosomes"]

    for sample_pair in sample_list:
        curr_samples = [file_name.strip() for file_name in sample_pair.get("filename").split(",")]
        fastq_end1, file_suffix, is_zipped = _separate_file_suffix(curr_samples[0])

        curr_output_address = output_address + "/{}/{}/{}".format(project_name, workflow, fastq_end1)

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
            input_address = curr_output_address

        #some precursors are never zipped
        if not can_be_zipped:
            is_zipped = "False"

        yield [project_name, workflow, file_suffix, ROOT_DIR, fastq_end1, fastq_end2, 
                input_address, curr_output_address, log_dir, is_zipped]

def _by_group_argument_generator(project_name, workflow, group_list, input_address, output_address, config_dictionary, log_dir):
    download_suffix = config_dictionary["download_suffix"]
    input_is_output = config_dictionary["input_is_output"]
    uses_chromosomes = config_dictionary["uses_chromosomes"]
    can_be_zipped = config_dictionary["can_be_zipped"]

    for group in group_list:
        samples = " ".join([x[0] for x in group_list[group]])
        file_suffix = group_list[group][0][1]
        is_zipped = group_list[group][0][2]
        curr_output_address = output_address + "/{}/{}/{}".format(project_name, workflow, group)

        if download_suffix:
            if uses_chromosomes:
                file_suffix = download_suffix
            else:
                file_suffix = download_suffix.format(file_suffix)

        if input_is_output:
            input_address = output_address + "/{}/{}".format(project_name, workflow)

        if not can_be_zipped:
            is_zipped = "False"

        yield [project_name, workflow, file_suffix, ROOT_DIR, group, "NULL", input_address,
                curr_output_address, log_dir, is_zipped, samples]

def _by_pair_argument_generator(project_name, workflow, group_list, pair_list, input_address, output_address, config_dictionary, log_dir):
    download_suffix = config_dictionary["download_suffix"]
    input_is_output = config_dictionary["input_is_output"]
    uses_chromosomes = config_dictionary["uses_chromosomes"]
    can_be_zipped = config_dictionary["can_be_zipped"]

    for group in group_list:
        first_sample, second_sample = [x[0] for x in group_list[group]]
        if pair_list.get(first_sample, None):
            normal_sample = first_sample
            tumor_sample = second_sample
        elif pair_list.get(second_sample, None):
            normal_sample = second_sample
            tumor_sample = first_sample 
        else:
            continue

        file_suffix = group_list[group][0][1]
        is_zipped = group_list[group][0][2]
        curr_output_address = output_address + "/{}/{}/{}".format(project_name, workflow, normal_sample)

        if download_suffix:
            if uses_chromosomes:
                file_suffix = download_suffix
            else:
                file_suffix = download_suffix.format(file_suffix)

        if not can_be_zipped:
            is_zipped = "False"

        if input_is_output:
            input_address = output_address + "/{}/{}".format(project_name, workflow)


        yield [project_name, workflow, file_suffix, ROOT_DIR, normal_sample, tumor_sample, input_address,
                curr_output_address, log_dir, is_zipped]



def _by_all_samples_argument_generator(project_name, workflow, sample_list, output_address, config_dictionary, log_dir):
    download_suffix = config_dictionary["download_suffix"]
    input_is_output = config_dictionary["input_is_output"]
    can_be_zipped = config_dictionary["can_be_zipped"]

    first_sample = sample_list[0].get("filename").split(",")[0].strip()
    fastq_end1, file_suffix, is_zipped = _separate_file_suffix(first_sample)

    # turn download suffix to an empty string when it's a Nonetype
    if download_suffix:
        if uses_chromosomes:
            file_suffix = download_suffix
        else:
            file_suffix - download_suffix.format(file_suffix)

    curr_output_address = output_address + "/{}/{}".format(project_name, workflow)

    if input_is_output:
        input_address = curr_output_address

    if not can_be_zipped:
        is_zipped = "False"

    samples = []

    for sample in sample_list:
        samples.append(sample.get("description"))

    samples = " ".join(samples)

    return [project_name, workflow, file_suffix, ROOT_DIR, "NULL", "NULL", input_address, 
            curr_output_address, log_dir, is_zipped, samples]


if __name__ == "__main__":
    run_analysis(sys.argv[1], sys.argv[2], sys.argv[3])
