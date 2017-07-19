__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import os
from util import YamlFileMaker
from util import QstatParser
from cfnCluster import ConnectionManager
import sys

workspace = "/shared/workspace/WGSPipeline/"
#log_dir = "/shared/workspace/data_archive/DNASeq/{}/logs"
log_dir = "/home/ec2-user/{}"

## executing WGS pipeline with the specific yaml file
def execute(ssh_client, project_name, analysis_steps, s3_input_files_address,
                   sample_list, group_name, s3_output_files_address, email):
    yaml_file = project_name + ".yaml"

    global log_dir
    log_dir = log_dir.format(project_name)

    print("making the yaml files ...")
    files = YamlFileMaker.make_wgs_yaml_files(yaml_file, project_name, analysis_steps, s3_input_files_address,
                   sample_list, group_name, s3_output_files_address, "hg19", "NA")

    print("copying yaml files to remote master node...")
    for f in files:
        ConnectionManager.copy_file(ssh_client, f, workspace + "yaml_examples")
        os.remove(f)

    files = [workspace + "yaml_examples/" + f for f in files]

    #if not email == "":

    print("executing pipeline...")
    ConnectionManager.execute_command(ssh_client, "qsub -e /dev/null -o /dev/null " + workspace + "scripts/" + "submit.sh "
             + ",".join(files) + " " + log_dir)


## checking your jobs status
def check_status(ssh_client, job_name):
    print("checking processing status")
    qstat = ConnectionManager.execute_command(ssh_client, "qstat")

    job_ids = QstatParser.get_job_ids(qstat)
    job_details = [ConnectionManager.execute_command(ssh_client, 
        "qstat -j %s" % x[0]) for x in job_ids]

    job_info = [job_ids[x] + [job_details[x]] for x in range(len(job_ids))]

    global log_dir
    logs = ConnectionManager.list_dir(ssh_client, log_dir)

    QstatParser.parse_qstat(job_info, job_name, logs)

## checking your jobs status
def check_jobs_status(ssh_client):
    print("checking jobs status")
    ConnectionManager.execute_command(ssh_client, "qstat")

## checking your host status
def check_host_status(ssh_client):
    print("checking qhost status")
    ConnectionManager.execute_command(ssh_client, "qhost")
