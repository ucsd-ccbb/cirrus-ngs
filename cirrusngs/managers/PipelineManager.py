__doc__="""
This module contains the main pipeline execution method on the 
local side. It is shared between all pipelines and workflows with
differing parameters.
"""

import os
from cirrusngs.util import YamlFileMaker
from cirrusngs.managers import ConnectionManager
from cirrusngs.util import QstatParser
import subprocess

workspace = "/shared/workspace/Pipelines/"

# executing a specified pipeline with the specific yaml file
def execute(pipeline, ssh_client, project_name, workflow, analysis_steps, s3_input_files_address,
        sample_list, group_list, s3_output_files_address, genome, style, pairs_list):
    """Executes a pipeline.

    The main local side function for executing a pipeline with all user inputs to jupyter notebook.
    Calls the run.sh shell script on the cluster head node using nohup after creating
    a yaml file summarizing user input and uploaded that file to the cluster.

    Args:
        pipeline: name of the pipeline to be run, supported pipelines can be found in CirrusAddons notebook
        ssh_client: a paramiko SSHClient object that connects to the cluster where analysis is run
        project_name: name of the current project, <project_name>.yaml contains all user input to notebook
        workflow: name of the workflow to be run, supported workflows can be found in CirrusAddons notebook
        analysis_steps: set of analysis steps to be run, supported steps can be found in pipeline's notebook
        s3_input_files_address: s3 bucket containing all fastq files for project
        sample_list: list of dictionaries with sample info for each sample
        group_list: list of all groups, shares indices with sample_list (sample_list[0] is in group[0], etc)
        s3_output_files_address: root s3 bucket where analysis results should be uploaded
        genome: reference genome to be used, supported genomes can be found in pipeline's notebook
        style: only for ChIPSeq homer workflow, can be "factor" or "histone"
        pairs_list: dictionary with keys=normal samples, values=experimental samples
            for ChIPSeq the keys=ChIP samples, values=corresponding input regularization samples

    Returns:
        None
    """
    yaml_file = project_name + ".yaml"

    if s3_output_files_address.endswith("/"):
        s3_output_files_address = s3_output_files_address[:-1]
    if s3_input_files_address.endswith("/"):
        s3_input_files_address = s3_input_files_address[:-1]

    logs_dir = "/shared/workspace/logs/{}/{}/{}".format(pipeline, workflow, project_name)

    print("making the yaml file...")
    YamlFileMaker.make_yaml_file(yaml_file, pipeline, project_name, workflow, analysis_steps, s3_input_files_address,
                                 sample_list, group_list, s3_output_files_address, genome, style, pairs_list)

    print("copying yaml file to remote master node...")

    # Make sure remote directory exists
    remote_dir = workspace + "yaml_files/" + pipeline + "/" + workflow
    ssh_client.exec_command("mkdir -p " + remote_dir)
    
    ConnectionManager.copy_file(ssh_client, yaml_file, "{}yaml_files/{}/{}".format(workspace, pipeline, workflow))

    # Remove the local yaml file
    os.remove(yaml_file)

    print("executing pipeline...")
    
    ConnectionManager.execute_command(ssh_client, "nohup bash " + workspace + "scripts/run.sh "
                                      + workspace + "yaml_files/{}/{}/{} ".format(pipeline, workflow, yaml_file)
                                      + logs_dir + " " + pipeline+"_"+workflow)

def check_status(ssh_client, step_name, pipeline, workflow, project_name,analysis_steps,verbose=False):
    QstatParser.check_status(ssh_client, step_name, pipeline, workflow, project_name,analysis_steps,verbose)

def stop_pipeline(ssh_client):
    ConnectionManager.execute_command(ssh_client, "bash /shared/workspace/Pipelines/util/stop_pipeline.sh")
