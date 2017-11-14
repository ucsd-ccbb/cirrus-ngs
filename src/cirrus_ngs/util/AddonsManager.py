from cfnCluster import ConnectionManager
from ast import literal_eval
from difflib import get_close_matches
import yaml

def get_scripts_dict(ssh_client):
    return literal_eval(ConnectionManager.execute_command(ssh_client, "python /shared/workspace/Pipelines/util/GetScripts.py"))

def get_all_pipeline_names(script_dict):
    return [pipeline for pipeline in script_dict if not pipeline == "All Pipelines"]

def get_workflows_in_pipeline(script_dict, pipeline):
    if pipeline == "all":
        print()
        return {curr_pipeline:[workflow for workflow in script_dict[curr_pipeline] if not workflow == "All Workflows"] 
                for curr_pipeline in script_dict if not curr_pipeline == "All Pipelines"}
    return [workflow for workflow in script_dict[pipeline] if not workflow == "All Workflows"]

def get_scripts(script_dict, pipeline, workflow):
    if pipeline == "all":
        return script_dict

    if workflow == "all":
        return {**script_dict[pipeline], **{"All Pipelines":script_dict["All Pipelines"]}}

    return {workflow:script_dict[pipeline][workflow], 
            "All {} Workflows".format(pipeline):script_dict[pipeline]["All Workflows"],
            "All Pipelines":script_dict["All Pipelines"]}

def get_steps_calling_script(ssh_client, script_dict, script_name):
    tools_conf = yaml.load(ConnectionManager.execute_command(ssh_client, "cat /shared/workspace/Pipelines/config/tools.yaml"))

    pipe_work_dict = {pipeline:get_workflows_in_pipeline(script_dict, pipeline) for pipeline in get_all_pipeline_names(script_dict)}



    result = "{} called from:\n".format(script_name)
    script_name = script_name.replace(".sh", "")
    for step in tools_conf:
        dirs = tools_conf[step]["script_path"].split("/")
        
        in_strings = ["in all Pipelines", "in all Workflows in {}", "in the {} {} workflow"]
        in_string = in_strings[len(dirs)-1].format(*dirs[:-1])

        if dirs[-1] == script_name:
            result += "{} {}\n".format(step, in_string)

    return result




def cat_script(ssh_client, script_dict, pipeline, workflow, script_name):
    if script_name in script_dict[pipeline][workflow]:
        return "Workflow Specific", ConnectionManager.execute_command(ssh_client, "cat /shared/workspace/Pipelines/scripts/{}/{}/{}".format(pipeline, workflow, script_name))
    elif script_name in script_dict[pipeline]["All Workflows"]:
        return "Pipeline Specific", ConnectionManager.execute_command(ssh_client, "cat /shared/workspace/Pipelines/scripts/{}/{}".format(pipeline, script_name))
    elif script_name in script_dict["All Pipelines"]:
        return "All Pipelines", ConnectionManager.execute_command(ssh_client, "cat /shared/workspace/Pipelines/scripts/{}".format(script_name))
    else:
        return "This script doesn't exist anywhere"


def get_software_dict(ssh_client):
    return literal_eval(ConnectionManager.execute_command(ssh_client, "python /shared/workspace/Pipelines/util/GetSoftware.py"))

def check_tool_is_installed(software_dict, tool):
    for key in software_dict:
        if key.lower() == tool.lower():
            return "{} installed with version(s):\n{}".format(key, "\n".join(software_dict[key]))

    recommendations = get_close_matches(tool.lower(), [k.lower() for k in software_dict.keys()])

    if recommendations:
        recommendations = " or ".join([k for k in software_dict.keys() if k.lower() in recommendations])
        return "Did you mean {}?".format(recommendations)
    else:
        return "{} not installed".format(tool)

