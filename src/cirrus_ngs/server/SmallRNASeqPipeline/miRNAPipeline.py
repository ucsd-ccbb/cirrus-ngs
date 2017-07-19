__author__ = 'Mengyi Liu<mel097@ucsd.edu>'

import sys
import subprocess
import PBSTracker
import YamlFileReader

workspace = "/shared/workspace/SmallRNASeqPipeline/"
scripts = "/shared/workspace/SmallRNASeqPipeline/scripts/"
data_dir = "/shared/workspace/data_archive/SmallRNASeq"  # for logs and download file


def runPipeline(yaml_file):
    
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
    
    
    # when only one sample
    if len(sample_list) < 2:
        download(project_name, file_suffix, data_dir, sample_list[0].get("filename")[:sample_list[0].get("filename").find(file_suffix)],
                 "NULL", s3_input_address)
    # more than one sample
    else:
        download(project_name, file_suffix, data_dir, sample_list[0].get("filename")[:sample_list[0].get("filename").find(file_suffix)],
                 sample_list[1].get("filename")[:sample_list[1].get("filename").find(file_suffix)], s3_input_address)

    # call runFastQC
    if "fastqc" in analysis_steps:
        if len(sample_list) < 2:
            runFastQC(project_name, file_suffix, data_dir,  sample_list[0].get("filename")[:sample_list[0].get("filename").find(file_suffix)],
                      "NULL", "hg19")
        else:
            runFastQC(project_name, file_suffix, data_dir,
                      sample_list[0].get("filename")[:sample_list[0].get("filename").find(file_suffix)],
                      sample_list[1].get("filename")[:sample_list[1].get("filename").find(file_suffix)], "hg19")


## download fastq files
def download(project_name, file_suffix, data_dir, file1_name, file2_name, s3_input_address):
    subprocess.call(['qsub', scripts + 'download.sh', project_name, file_suffix, data_dir, file1_name, file2_name, s3_input_address])
    
    PBSTracker.trackPBSQueue(1, "download")


## run fastqc
def runFastQC(project_name, file_suffix, data_dir, file1, file2, genome):
    
    print("executing fastqc...")
    print("Project name: " + project_name)
    
    subprocess.call(['qsub', scripts + "fastqc.sh",
                     project_name, file_suffix, data_dir, file1, file2, genome])
                     
    PBSTracker.trackPBSQueue(1, "fastqc")


if __name__ == "__main__":
    runPipeline(sys.argv[1])
