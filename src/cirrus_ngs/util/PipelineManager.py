import os
from util import YamlFileMaker
from cfnCluster import ConnectionManager
from util import QstatParser

workspace = "/shared/workspace/Pipelines/"


# executing a specified pipeline with the specific yaml file
def execute(pipeline, ssh_client, project_name, workflow, analysis_steps, s3_input_files_address,
        sample_list, group_list, s3_output_files_address, genome, style, pairs_list):
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
    ConnectionManager.copy_file(ssh_client, yaml_file, "{}yaml_files/{}/{}".format(workspace, pipeline, workflow))

    # Remove the local yaml file
    os.remove(yaml_file)

    print("executing pipeline...")
    
    ConnectionManager.execute_command(ssh_client,
                                      "nohup bash " + workspace + "scripts/run.sh "
                                      + workspace + "yaml_files/{}/{}/{} ".format(pipeline, workflow, yaml_file)
                                      + logs_dir + " " + pipeline+"_"+workflow)

def check_status(ssh_client, job_name):
    print("checking status of jobs...")
    qstat = ConnectionManager.execute_command(ssh_client, "qstat")
    logs = ConnectionManager.list_dir(ssh_client, logs_dir)
    QstatParser.parse_qstat(qstat, logs, job_name)
