import os
from util import YamlFileMaker
from cfnCluster import ConnectionManager

workspace = "/shared/workspace/Pipelines/"
logs_dir = "/shared/workspace/logs/{}/{}"


# executing a specified pipeline with the specific yaml file
def execute(pipeline, ssh_client, project_name, analysis_steps, s3_input_files_address,
            sample_list, group_list, s3_output_files_address, genome, style):
    yaml_file = project_name + ".yaml"
    # specify the log directory
    global logs_dir
    logs_dir = logs_dir.format(pipeline, project_name)

    print("making the yaml file...")
    YamlFileMaker.make_yaml_file(yaml_file, project_name, analysis_steps, s3_input_files_address,
                                 sample_list, group_list, s3_output_files_address, genome, style)

    print("copying yaml file to remote master node...")
    ConnectionManager.copy_file(ssh_client, yaml_file, workspace + "yaml_examples")

    # Remove the local yaml file
    os.remove(yaml_file)

    print("executing pipeline...")
    ConnectionManager.execute_command(ssh_client,
                                      "qsub -o /dev/null -e /dev/null " + workspace + "scripts/run.sh "
                                      + workspace + "yaml_examples/" + yaml_file + " " + logs_dir + " "
                                      + pipeline)
