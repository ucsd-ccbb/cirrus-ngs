__doc__ = """
This module contains a library of functions used in the
CirrusAddons jupyter notebook. It contains utility functions 
that can be used to check content on the cluster from jupyter.
"""
__author__ = "Mustafa Guler"

import os
from cfnCluster import ConnectionManager
from ast import literal_eval
from difflib import get_close_matches
import yaml
from pygments import highlight
from pygments.lexers import BashLexer 
from pygments.formatters import HtmlFormatter
from IPython.core.display import display, HTML

shell_script_template = """#!/bin/bash

project_name=$1
workflow=$2
file_suffix=$3
root_dir=$4
fastq_end1=$5
fastq_end2=$6
input_address=$7
output_address=$8
log_dir=$9
is_zipped=${10}
{EXTRA ARGUMENTS HERE}

#logging
log_dir=$log_dir/fastq_end1
mkdir -p $log_dir
log_file=$log_dir/{LOG FILE NAME HERE}
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$workflow/$fastq_end1
mkdir -p $workspace

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME $status_file

##DOWNLOAD##
if [ ! -f $workspace/$fastq_end1$file_suffix ]
then
    download_suffix=$file_suffix

    if [ "$is_zipped" == "True" ]
    then
        download_suffix=$file_suffix".gz"
    fi

    check_exit_status "aws s3 cp $input_address/$fastq_end1$download_suffix $workspace/" $JOB_NAME $status_file
    gunzip -q $workspace/$fastq_end1$download_suffix

    if [ "$fastq_end2" != "NULL" ]
    then
        check_exit_status "aws s3 cp $input_address/$fastq_end2$download_suffix $workspace/" $JOB_NAME $status_file
        gunzip -q $workspace/$fastq_end2$download_suffix
    fi
fi
##END_DOWNLOAD##


##TOOL##
check_exit_status "{TOOL CALL HERE}" $JOB_NAME $status_file
##END_TOOL##


##UPLOAD##
#can use multiple includes if needed 
aws s3 cp $workspace $output_address/ --exclude "*" --include "{GLOB FOR OUTPUT FILES HERE}" --recursive
##END_UPLOAD##
"""

def get_scripts_dict(ssh_client):
    """
    Returns a dictionary representing all scripts on the cluster. Uses cluster's GetScripts.py.
    The dictionary has supported Pipelines and "All Pipelines" as keys and a second dict as values.
        The inner dictionaries have supported workflows for that pipeline and "All Workflows" as keys.
        The values for each is a list of shell script names
    This dictionary is the core data structure behind the "Check Scripts on Cluster" section in the
    CirrusAddons notebook

    input:
        ssh_client: a paramiko SSHClient obj
    """
    return literal_eval(ConnectionManager.execute_command(ssh_client, "python /shared/workspace/Pipelines/util/GetScripts.py"))

def get_all_pipeline_names(scripts_dict):
    """
    Returns the names of the supported pipelines on the cluster. 
    Does so by returning all keys in scripts dictionary except for "All Pipelines"

    input: 
        scripts_dict: scripts dictionary described in get_scripts_dict function
    """
    return [pipeline for pipeline in scripts_dict if not pipeline == "All Pipelines"]

def get_workflows_in_pipeline(scripts_dict, pipeline):
    """
    IF pipeline == "all":
        return a dictionary with keys=pipeline names, values=list of workflows in that pipeline 
    ELSE:
        return a list of workflows in a given pipeline

    input: 
        scripts_dict: scripts dictionary described in get_scripts_dict function
        pipeline: name of a supported target pipeline or "all". 
    """
    if pipeline == "all":
        print()
        return {curr_pipeline:[workflow for workflow in scripts_dict[curr_pipeline] if not workflow == "All Workflows"] 
                for curr_pipeline in scripts_dict if not curr_pipeline == "All Pipelines"}
    return [workflow for workflow in scripts_dict[pipeline] if not workflow == "All Workflows"]

def get_scripts(scripts_dict, pipeline, workflow):
    """
    input:
        scripts_dict: scripts dictionary described in get_scripts_dict function
        pipeline: name of a supported target pipeline or "all". 
        workflow: name of a supported target workflow or "all"

    If the pipeline parameter is "all", returns the full scripts_dict
    If it isn't "all" but workflow is "all", returns dictionary with keys=workflows in pipeline, values=scripts in workflow
        It also contains the "All Pipelines" list from the original scripts_dict
    If both pipeline and workflow aren't "all" returns dictionary with keys=target workflow, All <pipeline> Workflows, All Pipelines
        and values=list of scripts under each key
    """
    if pipeline == "all":
        return scripts_dict

    if workflow == "all":
        return {**scripts_dict[pipeline], **{"All Pipelines":scripts_dict["All Pipelines"]}}

    return {workflow:scripts_dict[pipeline][workflow], 
            "All {} Workflows".format(pipeline):scripts_dict[pipeline]["All Workflows"],
            "All Pipelines":scripts_dict["All Pipelines"]}

def get_steps_calling_script(ssh_client, scripts_dict, script_name):
    """
    Gets which steps call the given shell script. Does so by parsing workflow config yamls on the cluster. 
    Returns str that summarizes that information

    input:
        ssh_client: a paramiko SSHClient obj
        scripts_dict: scripts dictionary described in get_scripts_dict function
        script_name: name of target script including .sh extension
    """
    tools_conf = yaml.load(ConnectionManager.execute_command(ssh_client, "cat /shared/workspace/Pipelines/config/tools.yaml"))

    #keys=pipline names, values=list of workflows in that pipeline
    #excludes "All Pipelines" key from original scripts_dict
    pipe_work_dict = {pipeline:get_workflows_in_pipeline(scripts_dict, pipeline) for pipeline in get_all_pipeline_names(scripts_dict)}

    configs = []
    for pipeline in pipe_work_dict:
        for workflow in pipe_work_dict[pipeline]:
            configs.append(yaml.load(ConnectionManager.execute_command(ssh_client, "cat /shared/workspace/Pipelines/config/{0}/{0}_{1}.yaml".format(pipeline, workflow))))

    result = "{} called from:\n".format(script_name)
    script_name = script_name.replace(".sh", "")

    for conf in configs:
        for step in conf:
            if step == "steps":
                continue

            dirs = conf[step]["script_path"].split("/")
            
            #if length of dirs 1, then the script path must just contain script name
            #   therefore, script in all Pipelines
            #if length of dirs 2, then script path must have pipeline and script name
            #   therefore, script in all workflows for a specific pipeline
            #if length of dirs 3, then script path must have pipeline,workflow,script name
            #   therefore, script in specific workflow for specific pipeline
            in_strings = ["in all Pipelines", "in all Workflows in {}", "in the {} {} workflow"]
            in_string = in_strings[len(dirs)-1].format(*dirs[:-1])

            if dirs[-1] == script_name:
                result += "{} {}\n".format(step, in_string)
                if len(dirs) != 3:
                    return result

    return result


def get_step_config_dicts(ssh_client, scripts_dict, step_name):
    tools_conf = yaml.load(ConnectionManager.execute_command(ssh_client, "cat /shared/workspace/Pipelines/config/tools.yaml"))
    step_tools_conf = tools_conf[step_name]

    specific_confs_dict = literal_eval(ConnectionManager.execute_command(ssh_client, "python /shared/workspace/Pipelines/util/GetAllSpecificConfs.py"))

    step_spec_conf = {}

    for pipeline in specific_confs_dict:
        for workflow in specific_confs_dict[pipeline]:
            curr_conf = yaml.load(specific_confs_dict[pipeline][workflow])
            if step_name in curr_conf:
                step_spec_conf.update({(pipeline, workflow):curr_conf[step_name]})

    return step_tools_conf, step_spec_conf

def get_step_config(ssh_client, scripts_dict, step_name, step_tools_conf, step_spec_conf):
    step_conf = yaml.dump(step_tools_conf, default_flow_style=False)
    result = "\ntools.yaml configuration entry for {} step:\n{}\n\n".format(step_name, step_conf)

    for key in step_spec_conf:
        curr_entry = yaml.dump({step_name:step_spec_conf[key]}, default_flow_style=False)
        pipeline, workflow = key
        _,args = cat_script(ssh_client, scripts_dict, pipeline, workflow, step_tools_conf["script_path"].split("/")[-1] + ".sh")
        args = list(map(lambda x : x.split("=")[0], args.split("\n\n")[1].splitlines()[10:]))[:len(curr_entry.splitlines())-1]

        result += "{}_{}.yaml configuration entry for {} step:\n{}".format(pipeline, workflow, step_name, curr_entry)

        for ind,arg in enumerate(args):
            result += "Argument {} is {}, ".format(ind+1, arg)

        result = result.rstrip(", ")
        result += "\n\n"

    return result

def cat_script(ssh_client, scripts_dict, pipeline, workflow, script_name):
    """
    Returns tuple:
        (specificity of script, string representation of given script contents )
    Specificty in ["Workflow Specific", "Pipeline Specific", "All Pipelines"]
    Must specify which workflow and pipeline the script is in.

    input:
        ssh_client: a paramiko SSHClient obj
        scripts_dict: scripts dictionary described in get_scripts_dict function
        pipeline: name of a supported target pipeline
        workflow: name of a supported target workflow
        script_name: name of target script including .sh extension
    """
    if script_name in scripts_dict[pipeline][workflow]:
        return "Workflow Specific", ConnectionManager.execute_command(ssh_client, "cat /shared/workspace/Pipelines/scripts/{}/{}/{}".format(pipeline, workflow, script_name))
    elif script_name in scripts_dict[pipeline]["All Workflows"]:
        return "Pipeline Specific", ConnectionManager.execute_command(ssh_client, "cat /shared/workspace/Pipelines/scripts/{}/{}".format(pipeline, script_name))
    elif script_name in scripts_dict["All Pipelines"]:
        return "All Pipelines", ConnectionManager.execute_command(ssh_client, "cat /shared/workspace/Pipelines/scripts/{}".format(script_name))
    else:
        return "This script isn't called in the specified Pipeline/Workflow", ""

def show_script(str_script):
    html_vers = highlight(str_script, BashLexer(), HtmlFormatter(full=True))
    display(HTML(html_vers))


def edit_step_tools_config(ssh_client, new_step_tools_conf, step_name):
    tools_conf = yaml.load(ConnectionManager.execute_command(ssh_client, "cat /shared/workspace/Pipelines/config/tools.yaml"))
    tools_conf[step_name] = new_step_tools_conf
    with open("tools.yaml", "w+") as f:
        f.write(yaml.dump(tools_conf, default_flow_style=False))
    ConnectionManager.execute_command(ssh_client, "mv -n /shared/workspace/Pipelines/config/tools.yaml /shared/workspace/Pipelines/config/tools.yaml.BACKUP")
    ConnectionManager.copy_file(ssh_client, "{}/tools.yaml".format(os.getcwd()), "/shared/workspace/Pipelines/config/tools.yaml")

def edit_step_specific_config(ssh_client, pipeline, workflow, new_extra_bash_args, step_name):
    conf_file_name = "{}_{}.yaml".format(pipeline, workflow)
    spec_conf = yaml.load(ConnectionManager.execute_command(ssh_client, "cat /shared/workspace/Pipelines/config/{}/{}".format(pipeline, conf_file_name)))
    spec_conf[step_name] = new_extra_bash_args

    with open(conf_file_name, "w+") as f:
        f.write(yaml.dump(spec_conf, default_flow_style=False))
    ConnectionManager.execute_command(ssh_client, "mv -n /shared/workspace/Pipelines/config/{0}/{1} /shared/workspace/Pipelines/config/{0}/{1}.BACKUP".format(pipeline, conf_file_name))
    ConnectionManager.copy_file(ssh_client, "{}/{}".format(os.getcwd(), conf_file_name), "/shared/workspace/Pipelines/config/{}/{}".format(pipeline, conf_file_name))

def edit_script(ssh_client, scripts_dict, pipeline, workflow, script_name):
    _,script_text = cat_script(ssh_client, scripts_dict, pipeline, workflow, script_name)
    return "%%writefile {}\n{}".format(script_name, script_text)

def upload_script(ssh_client, pipeline, workflow, script_name):
    script_path_cluster = "/shared/workspace/Pipelines/scripts/"

    if pipeline == "all":
        script_path_cluster += script_name
    elif workflow == "all":
        script_path_cluster += "{}/{}".format(pipeline, script_name)
    else:
        script_path_cluster += "{}/{}/{}".format(pipeline, workflow, script_name)

    ConnectionManager.execute_command(ssh_client, "mv -n {0} {0}.BACKUP".format(script_path_cluster))
    ConnectionManager.copy_file(ssh_client, "{}/{}".format(os.getcwd(), script_name), script_path_cluster)

def restore_backups(ssh_client):
    ConnectionManager.execute_command(ssh_client, "python /shared/workspace/Pipelines/util/RestoreBackups.py")


def get_software_dict(ssh_client):
    """
    Returns a dictionary representing all software on the cluster. Uses cluster's GetSoftware.py.
    The dictionary has installed software names as keys and a list of the installed versions as values
    This dictionary is the core data structure behind the "Check Software on Cluster" section in the
    CirrusAddons notebook

    input:
        ssh_client: a paramiko SSHClient obj
    """
    return literal_eval(ConnectionManager.execute_command(ssh_client, "python /shared/workspace/Pipelines/util/GetSoftware.py"))

def check_tool_is_installed(software_dict, tool):
    """
    Returns str with tool's versions in it if the tool is installed on cluster, otherwise recommends other tools names if possible.
    If no recommendations availabe, returns str saying tools not installed

    input:
        software_dict: software dictionary described in get_software_dict function
        tool: name of a target tool/software to check installation of on cluster
    """
    for key in software_dict:
        if key.lower() == tool.lower():
            return "{} installed with version(s):\n{}".format(key, "\n".join(software_dict[key]))

    recommendations = get_close_matches(tool.lower(), [k.lower() for k in software_dict.keys()])

    if recommendations:
        recommendations = " or ".join([k for k in software_dict.keys() if k.lower() in recommendations])
        return "Did you mean {}?".format(recommendations)
    else:
        return "{} not installed".format(tool)

