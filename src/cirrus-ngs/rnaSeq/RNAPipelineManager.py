__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import os
import YamlFileMaker
from cfnCluster import ConnectionManager
from util import DesignFileLoader

workspace = "/shared/workspace/RNASeqPipeline"
data_dir = "/shared/workspace/data_archive/RNASeq"

## run all analysis from download, alignment, counting and differential calculation.
def run_analysis(ssh_client, workflow, project_name, analysis_steps,
                 s3_input_files_address, sample_list, group_list, s3_output_files_address):

    yaml_file = project_name + ".yaml"

    print "making the yaml file..."
    YamlFileMaker.make_yaml_file(yaml_file, workflow, project_name, analysis_steps, s3_input_files_address,
                   sample_list, group_list, s3_output_files_address)

    print "copying yaml file to remote master node..."
    ConnectionManager.copy_file(ssh_client, yaml_file, workspace + "/" + workflow + "/yaml_examples")

    ## Remove the local yaml file
    os.remove(yaml_file)

    print "executing pipeline..."
    ConnectionManager.execute_command(ssh_client, "sh " + workspace + "/run.sh "
                                      + workspace + "/" + workflow + "/yaml_examples/" + yaml_file)

## checking your jobs status
def check_processing_status(ssh_client):
    print "checking processing status"
    ConnectionManager.execute_command(ssh_client, "cat " + workspace + "/nohup.out")

## checking your jobs status
def check_jobs_status(ssh_client):
    print "checking jobs status"
    ConnectionManager.execute_command(ssh_client, "qstat")

## checking your host status
def check_host_status(ssh_client):
    print "checking qhost status"
    ConnectionManager.execute_command(ssh_client, "qhost")

if __name__ == '__main__':
    workflow = "star_htseq_workflow"
    analysis_steps = ["fastqc"]

    s3_input_files_address = "s3://ccbb-analysis/Guorong/jupyter-genomics/data_archive/test_data/ChiPSeq/RA2284"
    s3_output_files_address = "s3://ccbb-analysis/Guorong/jupyter-genomics/data_archive/analysis_data/ChiPSeq"
    project_name = "Sample_cDNA"
    ssh_client = ""
    design_file = "/Users/guorongxu/Desktop/workspace/projects/jupyter-genomics_bitbucket/data/awsCluster/rnaSeq_design_example.txt"

    sample_list, group_list = DesignFileLoader.load_design_file(design_file)

    run_analysis(ssh_client, workflow, project_name, analysis_steps,
                 s3_input_files_address, sample_list, group_list, s3_output_files_address)