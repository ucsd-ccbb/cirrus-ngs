#__doc__ =  """
#This utility module is used to create a yaml file summarizing the data
#a user put into a given pipeline's jupyter notebook. That yaml file
#is then moved to the cluster where analysis will be run
#"""
__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import time

def make_yaml_file(yaml_file, pipeline, project_name, workflow, analysis_steps, s3_input_files_address,
        sample_list, group_list, s3_output_files_address, genome, style, pairs_list):
    """Writes user input to a yaml file.

    The user information is taken from a pipeline's jupyter notebook. Yaml
    file is later used in the PipelineManager execute function.

    Args:
        yaml_file: A string name of file being written to
        pipeline: A string name of the pipeline being run
        project_name: A string name of the project being run
        workflow: A string name of the workflow being run under given pipeline
        analysis_steps: A set of strings corresponding to possible analysis steps 
            for the workflow. See notebooks for possible steps.
        s3_input_files_address: A string user-given s3 bucket path to 
            input files for analysis
        sample_list: from DesignFileLoader, see load_design_file function output
        group_list: from DesignFileLoader, see load_design_file function output
        s3_output_files_address: A string user-given s3 bucket path to 
            upload analysis output 
        genome: A string name of reference genome being used. 
            See pipeline notebook for allowed values
        style: Only for ChIPSeq, a string that is either "factor" or "histone"
        pairs_list: from DesignFileLoader, see load_design_file function output

    Returns:
        None
    """
    filewriter = open(yaml_file, "w")
    filewriter.write("project: " + project_name + "\n")
    filewriter.write("pipeline: " + pipeline + "\n")
    filewriter.write("workflow: " + workflow + "\n")

    filewriter.write("analysis:\n")
    for analysis in analysis_steps:
        filewriter.write("  - {}\n".format(analysis))

    filewriter.write("date: " + time.strftime("%Y/%m/%d") + "\n")
    filewriter.write("upload: " + s3_output_files_address + "\n")
    filewriter.write("download: " + s3_input_files_address + "\n")

    filewriter.write("genome: " + genome + "\n")
    filewriter.write("style: " + style + "\n")
    filewriter.write("pairs:\n")
    if not pairs_list:
        filewriter.write("  {}\n")
    else:
        for normal, tumor in pairs_list.items():
            filewriter.write("  {}: {}\n".format(normal, tumor))

    
    filewriter.write("sample:\n")
    for index, sample in enumerate(sample_list):
        filename_string = "- filename: "

        if len(sample) == 1:
            filewriter.write(filename_string + sample[0].rstrip() + "\n")
            filewriter.write("  description: " + sample[0].split(".")[0] + "\n")
            filewriter.write("  group: " + group_list[index] + "\n")
        elif len(sample) == 2:
            for sample_name in sample:
                filename_string = filename_string + sample_name + ", "
            filewriter.write(filename_string[:-2] + "\n")
            filewriter.write("  description: " + sample[0].split(".")[0] + "\n")
            filewriter.write("  group: " + group_list[index] + "\n")

    filewriter.close()
