__author__ = 'Mengyi Liu<mel097@ucsd.edu>'
# this file is in: /shared/workspace/SmallRNASeqPipeline

import sys
import subprocess
import PBSTracker
import YamlFileReader

workspace = "/shared/workspace/SmallRNASeqPipeline/"
scripts = "/shared/workspace/SmallRNASeqPipeline/scripts/"
data_dir = "/shared/workspace/data_archive/SmallRNASeq"  # for logs and download file


def run_analysis(yaml_file):

    documents = YamlFileReader.parse_yaml_file(yaml_file)

    project_name = documents.get("project")
    analysis_steps = documents.get("analysis")
    s3_input_address = documents.get("upload")
    sample_list = documents.get("sample")

    # print "sample_list[0]: ", sample_list[0]  # this is a dictionary
    # print "Length of sample list: ", len(sample_list)
    # print "filename: ", sample_list[0].get("filename")

    # set file suffix
    file_suffix = ".fastq"
    if sample_list[0].get("filename").find(".fq") > -1:
        file_suffix = ".fq"

    # call download
    # when only one sample
    if len(sample_list) < 2:
        download(project_name, file_suffix, data_dir, sample_list[0].get("filename")[:sample_list[0]
                 .get("filename").find(file_suffix)], "NULL", s3_input_address)
    # more than one sample
    else:
        download(project_name, file_suffix, data_dir, sample_list[0].get("filename")[:sample_list[0]
                 .get("filename").find(file_suffix)],
                 sample_list[1].get("filename")[:sample_list[1].get("filename").find(file_suffix)], s3_input_address)

    # call run_fastqc
    if "fastqc" in analysis_steps:
        if len(sample_list) < 2:
            run_fastqc(project_name, file_suffix, data_dir, sample_list[0].get("filename")[:sample_list[0]
                       .get("filename").find(file_suffix)], "NULL", "hg19")
        else:
            run_fastqc(project_name, file_suffix, data_dir,
                       sample_list[0].get("filename")[:sample_list[0].get("filename").find(file_suffix)],
                       sample_list[1].get("filename")[:sample_list[1].get("filename").find(file_suffix)], "hg19")

    # call run_trimmomatic
    if "trim" in analysis_steps:
        if len(sample_list) < 2:
            run_trimmomatic(project_name, file_suffix, data_dir,
                            sample_list[0].get("filename")[:sample_list[0].get("filename").find(file_suffix)], "NULL")
        else:
            run_trimmomatic(project_name, file_suffix, data_dir,
                            sample_list[0].get("filename")[:sample_list[0].get("filename").find(file_suffix)],
                            sample_list[1].get("filename")[:sample_list[1].get("filename").find(file_suffix)])


# download fastq files
# file1_name and file2_name are the sample(s), in case there's only one sample, file2_name is NULL, same below
def download(project_name, file_suffix, file_dir, file1_name, file2_name, s3_input_address):
    print("downloading...")
    subprocess.call(['qsub', scripts + 'download.sh', project_name, file_suffix, file_dir, file1_name, file2_name,
                     s3_input_address])

    PBSTracker.trackPBSQueue(1, "download")


# run fastqc
def run_fastqc(project_name, file_suffix, file_dir, file1_name, file2_name, genome):
    print("executing fastqc...")
    # print("Project name: " + project_name)
    subprocess.call(['qsub', scripts + "fastqc.sh",
                    project_name, file_suffix, file_dir, file1_name, file2_name, genome])

    PBSTracker.trackPBSQueue(1, "fastqc")


# run trimmomatic
def run_trimmomatic(project_name, file_suffix, file_dir, file1_name, file2_name):
    print("executing trimmomatic...")
    print("file1: ", file1_name)
    print("file2: ", file2_name)
    subprocess.call(['qsub', scripts + "trim.sh",
                    project_name, file_suffix, file_dir, file1_name, file2_name])

    PBSTracker.trackPBSQueue(0.25, "trim")


if __name__ == "__main__":
    run_analysis(sys.argv[1])
