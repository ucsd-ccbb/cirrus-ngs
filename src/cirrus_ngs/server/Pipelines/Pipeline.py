import sys
import os
import subprocess
from util import PBSTracker
import re
import yaml

ROOT_DIR = "/scratch"
SCRIPTS = "/shared/workspace/Pipelines/scripts/"

##run the actual pipeline
def run_analysis(yaml_file, log_dir, pipeline_config_file):
    """Runs a full pipeline.
    
    The function called from this module's main, which is 
    subsequently called from run.sh. User information comes
    from a yaml file uploaded by the PipelineManager. Additional
    configuration comes from the configuration files.

    Args:
        yaml_file: string absolute path to the yaml file uploaded by PipelineManager
        log_dir: string absolute path to directory in which logs should be stored
        pipeline_config_file: string absolute path to the specific configuration
            file for this pipeline and workflow.

    Returns:
        None
    """
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

    # set environment variables
    os.environ["style"] = documents.get("style")
    genome = documents.get("genome")
    os.environ["genome"] = genome

    #for chipseq
    if os.environ["style"] == "histone":
        os.environ["style_ext"] = ".regions.txt"
    elif os.environ["style"] == "factor":
        os.environ["style_ext"] = ".peaks.txt"

    
    os.environ["genome_fasta"] = os.environ[genome + "_fasta"]
    os.environ["genome_fai"] = os.environ.get(genome + "_fai", "")
    os.environ["dbsnp"] = os.environ.get(genome + "_dbsnp", "")
    os.environ["bwa_index"] = os.environ.get(genome + "_bwa_index", "")
    os.environ["bowtie_index"] = os.environ.get(genome + "_bowtie_index", "")
    os.environ["mills"] = os.environ.get(genome + "_mills", "")
    os.environ["hapmap"] = os.environ.get(genome + "_hapmap", "")
    os.environ["omni"] = os.environ.get(genome + "_omni", "")
    os.environ["snps_1000G"] = os.environ.get(genome + "_snps_1000G", "")
    os.environ["indels_1000G"] = os.environ.get(genome + "_indels_1000G" , "")
    os.environ["indels"] = os.environ.get(genome + "_indels", "")
    os.environ["cosmic"] = os.environ.get(genome + "_cosmic", "")
    os.environ["chromosome_list"] = os.environ.get(genome + "_chromosome_list", "")
    os.environ["genome_gtf"] = os.environ.get(genome + "_gtf", "")
    os.environ["STAR_index"] = os.environ.get(genome + "_STAR_index", "")
    os.environ["genome_bowtie2_index"] = os.environ.get(genome + "_bowtie2_index", "")


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

def run_tool(tool_config_dict, extra_bash_args, project_name, workflow, sample_list, input_address, output_address, group_list, pair_list, log_dir):
    """Runs a single step in the current pipeline.

    This function can run any step in any pipeline/workflow based
    on the configuration entries in tools.yaml and the pipeline
    specific configuration file for that step. Generates arguments
    for each sample with a different generator depending
    on how the step's tools.yaml entry. Eventually submits
    the shell script for the current step as a job for
    each sample.

    Args:
        tool_config_dict: configuration dictionary for this tool, parsed from tools.yaml
        extra_bash_args: additional arguments list to the shell script for this tool, taken from pipeline-specific configuration
        project_name: name of the project this tool is running in
        workflow: the workflow that this project is using
        sample_list: a list of dictionaries taken from the project_name.yaml file
        input_address: the s3 input address, taken from the project_name.yaml file
        output_address: the s3 input address, taken from the project_name.yaml file
        group_list: dictionary with key=group and value=list of tuples from separate_file_suffix
        pair_list: dictionary with key=normal_sample and value=associated tumor sample
        log_dir: the path to the log directory this tool

    Returns:
        None
    """
    #num_threads used to request nodes with a specific amount of threads for step
    #first extra argument should always be number of threads
    if len(extra_bash_args) > 0:
        num_threads = extra_bash_args[0]
    else:
        num_threads = 1

    #chromosome lists are determined by the genome name from the software.conf configuration file
    CHROMOSOME_LIST = os.environ["chromosome_list"].split()

    #contains all qsub flags and the name of the shell script for this step
    subprocess_call_list = ["qsub", "-V", "-o", "/dev/null", "-e", "/dev/null", "-pe", "smp", str(num_threads),
            SCRIPTS + tool_config_dict["script_path"] + ".sh"]
    
    #extra arguments to the current shell script outside of the 10 necessary ones (see shell script format documentation)
    #if not empty, first extra argument must be a number representing number of threads used by script
    extra_bash_args = list(map(str, extra_bash_args))

    #runs tool on every sample in project in one step 
    if tool_config_dict.get("all_samples", False):
        all_sample_arguments = _by_all_samples_argument_generator(project_name, 
                workflow, sample_list, input_address, output_address, tool_config_dict, log_dir)
        if tool_config_dict["uses_chromosomes"]:
            original_suffix = all_sample_arguments[2]
            for chromosome in CHROMOSOME_LIST:
                all_sample_arguments[2] = original_suffix.format(chromosome)
                subprocess.call(subprocess_call_list + all_sample_arguments + extra_bash_args)
        else:
            subprocess.call(subprocess_call_list + all_sample_arguments + extra_bash_args)
        return

    #run tool on pairs in project
    if tool_config_dict.get("by_pair", False) and os.environ["pairs_exist"] == "True":
        for pair_arguments in _by_pair_argument_generator(project_name, workflow, 
                group_list, pair_list, input_address, output_address, tool_config_dict, log_dir):
            if tool_config_dict["uses_chromosomes"]:
                original_suffix = pair_arguments[2]
                for chromosome in CHROMOSOME_LIST:
                    pair_arguments[2] = original_suffix.format(chromosome)
                    subprocess.call(subprocess_call_list + pair_arguments + extra_bash_args + [chromosome])
            else:
                subprocess.call(subprocess_call_list + pair_arguments + extra_bash_args)
        
        PBSTracker.trackPBSQueue(1, tool_config_dict["script_path"])
        return

    #runs tool on samples in each group
    if tool_config_dict.get("by_group", False):
        for group_arguments in _by_group_argument_generator(project_name, workflow, group_list, input_address, output_address, tool_config_dict, log_dir):
            if tool_config_dict["uses_chromosomes"]:
                original_suffix = group_arguments[2]
                for chromosome in CHROMOSOME_LIST:
                    group_arguments[2] = original_suffix.format(chromosome)
                    subprocess.call(subprocess_call_list + group_arguments + 
                            extra_bash_args + [chromosome])
            else:
                subprocess.call(subprocess_call_list + group_arguments + extra_bash_args)

        PBSTracker.trackPBSQueue(1, tool_config_dict["script_path"])
        return

    #run tool separately for each sample
    for curr_sample_arguments in _sample_argument_generator(project_name, workflow, sample_list, input_address, output_address, tool_config_dict, log_dir):
        #for tools that run on each chromosome the file suffix has the current chrom number added to it
        if tool_config_dict["uses_chromosomes"]:
            original_suffix = curr_sample_arguments[2]
            for chromosome in CHROMOSOME_LIST:
                curr_sample_arguments[2] = original_suffix.format(chromosome)
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
    """Extracts information about file extensions from fq/fastq files.

    FASTQ files can come with .fq or .fastq extensions, this
    function automatically extracts those extensions. Additionally, 
    this function determines if the fastq file was gzipped.

    Args:
        sample_file: a string name of a file, not the path

    Returns:
        a tuple of form:
        (
            string file name without extensions, 
            the .fq or .fastq extension,
            string "True" or "False" representing if file was gzipped, 
        )
    """
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
    """Generator that yields list of shell script arguments for sample-by-sample tool runs.

    For each sample in the project generates the shells script
    arguments and yields them. These arguments will be passed to the
    shell script when it's submitted.

    Args:
        project_name: string name of the project being run
        workflow: string workflow this run belongs to
        sample_list: a list of dictionaries taken from the project_name.yaml file
        input_address: the s3 input address, taken from the project_name.yaml file
        output_address: the s3 input address, taken from the project_name.yaml file
        config_dictionary: configuration dictionary for this tool, taken from tools.yaml
        log_dir: the path to the log directory this tool

    Yields:
        list
        [
            project name, 
            workflow, 
            file suffix used for downloading file, 
            the directory where outputs should be stored while running, 
            the forward reads file name (no ext), 
            the reverse reads file name (no ext, NULL if no reverse), 
            the s3 input address for precursor files,
            the s3 output address for uploading results, 
            the directory on the cluster for log files, 
            a str(bool) representing if the input files are gzipped or not
        ]
        for each sample
    """
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
    """Generator that yields list of shell script arguments for group-by-group tool runs.

    For each group in the project generates the shell script
    arguments and yields them. These arguments will be passed to the
    shell script when it's submitted.

    Args:
        project_name: string name of the project being run
        workflow: string workflow this run belongs to
        group_list: a dictionary with key=group name, value=list tuples for samples in that group from separate_file_suffix
        input_address: the s3 input address, taken from the project_name.yaml file
        output_address: the s3 input address, taken from the project_name.yaml file
        config_dictionary: configuration dictionary for this tool, taken from tools.yaml
        log_dir: the path to the log directory this tool

    Yields:
        list
        [
            project name, 
            workflow, 
            file suffix used for downloading file, 
            the directory where outputs should be stored while running, 
            the name of the current group,
            "NULL",
            the s3 input address for precursor files,
            the s3 output address for uploading results, 
            the directory on the cluster for log files, 
            a str(bool) representing if the input files are gzipped or not,
            a space delimited string with all samples in this group in it
        ]
        for each group in the project
    """
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
    """Generator that yields list of shell script arguments for pair-by-pair tool runs.

    For each pair in the project generates the shell script
    arguments and yields them. These arguments will be passed to the
    shell script when it's submitted. Pairs can be any two
    samples that are somehow related. Used in the bwa_mutect
    workflow in DNASeq for normal_vs_tumor comparisons as well
    as in ChIPSeq for input_vs_chip comparisons.

    Args:
        project_name: name of the project being run
        workflow: workflow this run belongs to
        group_list: a dictionary with key=group name, value=list tuples for samples in that group from separate_file_suffix
        pair_list: a dictionary with key=normal sample, value=tumor sample
        input_address: the s3 input address, taken from the project_name.yaml file
        output_address: the s3 input address, taken from the project_name.yaml file
        config_dictionary: configuration dictionary for this tool, taken from tools.yaml
        log_dir: the path to the log directory this tool

    Yields:
        list
        [
            project name, 
            workflow, 
            file suffix used for downloading file, 
            the directory where outputs should be stored while running, 
            the name of the normal sample,
            the name of the tumor sample,
            the s3 input address for precursor files,
            the s3 output address for uploading results, 
            the directory on the cluster for log files, 
            a str(bool) representing if the input files are gzipped or not,
        ]
        for each pair in the project
    """
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


def _by_all_samples_argument_generator(project_name, workflow, sample_list, input_address, output_address, config_dictionary, log_dir):
    """Generator that yields list of shell script arguments for tools that run on all samples at once.

    Generates one argument list that contains information 
    about all the samples in the project. These arguments will be passed 
    to the shell script when it's submitted. All sample runs can be
    useful for comparisons (multiqc, counting, etc)

    Args:
        project_name: name of the project being run
        workflow: workflow this run belongs to
        sample_list: a list of dictionaries taken from the project_name.yaml file
        input_address: the s3 input address, taken from the project_name.yaml file
        output_address: the s3 input address, taken from the project_name.yaml file
        config_dictionary: configuration dictionary for this tool, taken from tools.yaml
        log_dir: the path to the log directory this tool

    Yields:
        list
        [
            project name, 
            workflow, 
            file suffix used for downloading file, 
            the directory where outputs should be stored while running, 
            NULL, 
            NULL, 
            the s3 input address for precursor files,
            the s3 output address for uploading results, 
            the directory on the cluster for log files, 
            a str(bool) representing if the input files are gzipped or not,
            a space delimited string of all samples in this project
        ]
        exactly once
    """
    download_suffix = config_dictionary["download_suffix"]
    input_is_output = config_dictionary["input_is_output"]
    can_be_zipped = config_dictionary["can_be_zipped"]
    uses_chromosomes = config_dictionary["uses_chromosomes"]

    first_sample = sample_list[0].get("filename").split(",")[0].strip()
    fastq_end1, file_suffix, is_zipped = _separate_file_suffix(first_sample)

    # turn download suffix to an empty string when it's a Nonetype
    if download_suffix:
        if uses_chromosomes:
            file_suffix = download_suffix
        else:
            file_suffix = download_suffix.format(file_suffix)

    curr_output_address = output_address + "/{}/{}".format(project_name, workflow)

    if input_is_output:
        input_address = curr_output_address

    if not can_be_zipped:
        is_zipped = "False"

    samples = []

    for sample in sample_list:
        samples.append(sample.get("description"))

    samples = " ".join(samples)

    yield [project_name, workflow, file_suffix, ROOT_DIR, "NULL", "NULL", input_address, 
            curr_output_address, log_dir, is_zipped, samples]


if __name__ == "__main__":
    run_analysis(sys.argv[1], sys.argv[2], sys.argv[3])
