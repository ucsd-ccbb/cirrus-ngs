__doc__ = """
This module contains a library of functions used in the
CirrusAddons jupyter notebook. It contains utility functions 
that can be used to check content on the cluster from jupyter.
"""
__author__ = "Mustafa Guler"

from cfnCluster import ConnectionManager
from ast import literal_eval
from difflib import get_close_matches
import yaml

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
    Gets which steps call the given shell script. Does so by parsing tools.yaml on the cluster. 
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



    result = "{} called from:\n".format(script_name)
    script_name = script_name.replace(".sh", "")
    for step in tools_conf:
        dirs = tools_conf[step]["script_path"].split("/")
        
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

