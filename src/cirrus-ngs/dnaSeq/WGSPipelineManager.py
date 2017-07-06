__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import os
from util import YamlFileMaker
from util import QstatParser
from cfnCluster import ConnectionManager

workspace = "/shared/workspace/WGSPipeline/"

## executing WGS pipeline with the specific yaml file
def execute(ssh_client, project_name, analysis_steps, s3_input_files_address,
                   sample_list, group_name, s3_output_files_address, email):
    yaml_file = project_name + ".yaml"

    print("making the yaml file...")
    YamlFileMaker.make_yaml_file(yaml_file, project_name, analysis_steps, s3_input_files_address,
                   sample_list, group_name, s3_output_files_address, "hg19", "NA")

    print("copying yaml file to remote master node...")
    ConnectionManager.copy_file(ssh_client, yaml_file, workspace + "yaml_examples")

    ## Remove the local yaml file
    os.remove(yaml_file)

    #if not email == "":

    print("executing pipeline...")
    ConnectionManager.execute_command(ssh_client, "qsub " + workspace + "run.sh "
                                      + workspace + "yaml_examples/" + yaml_file)

## checking your jobs status
def check_processing_status(ssh_client, job_name):
    print("checking processing status")
    qstat = ConnectionManager.execute_command(ssh_client, "qstat")
    QstatParser.parse_qstat(qstat, job_name)


## checking your jobs status
def check_jobs_status(ssh_client):
    print("checking jobs status")
    ConnectionManager.execute_command(ssh_client, "qstat")

## checking your host status
def check_host_status(ssh_client):
    print("checking qhost status")
    ConnectionManager.execute_command(ssh_client, "qhost")
