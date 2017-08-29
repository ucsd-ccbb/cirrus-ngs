import os
from util import YamlFileMaker
from cfnCluster import ConnectionManager
from util import QstatParser

workspace = "/shared/workspace/Pipelines/"
logs_dir = "/shared/workspace/logs/{}/{}"


# executing a specified pipeline with the specific yaml file
def execute(pipeline, ssh_client, project_name, analysis_steps, s3_input_files_address,
            sample_list, group_list, s3_output_files_address, genome, style, mutect_pairs):
    yaml_file = project_name + ".yaml"

    # specify the log directory
    global logs_dir
    # for RNA-seq, needs directory for the specific workflow
    if pipeline.startswith("RNA"):
        general_pipeline = "RNASeq"
        workflow = pipeline[7:]
        logs_dir = logs_dir.format(general_pipeline, project_name + "/" + workflow)
    else:
        logs_dir = logs_dir.format(pipeline, project_name)

    print("making the yaml file...")
    YamlFileMaker.make_yaml_file(yaml_file, project_name, analysis_steps, s3_input_files_address,
                                 sample_list, group_list, s3_output_files_address, genome, style, mutect_pairs)


    print("copying yaml file to remote master node...")
    ConnectionManager.copy_file(ssh_client, yaml_file, workspace + "yaml_examples")

    # Remove the local yaml file
    os.remove(yaml_file)

    print("executing pipeline...")
    ConnectionManager.execute_command(ssh_client,
                                      "qsub -V -o /dev/null -e /dev/null " + workspace + "scripts/run.sh "
                                      + workspace + "yaml_examples/" + yaml_file + " " + logs_dir + " "
                                      + pipeline)


def check_status(ssh_client, job_name):
    print("checking status of jobs...")
    qstat = ConnectionManager.execute_command(ssh_client, "qstat")
    logs = ConnectionManager.list_dir(ssh_client, logs_dir)
    QstatParser.parse_qstat(qstat, logs, job_name)