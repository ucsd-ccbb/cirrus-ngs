import sys
import os
import subprocess
import PBSTracker
import YamlFileReader

root_dir = "/shared/workspace/WGSPipeline"
data_dir = "/shared/workspace/data_archive/DNASeq"
data_files_dir = ""

##run the actual pipeline
def run_analysis(yaml_file):
    documents = YamlFileReader.parse_yaml_file(yaml_file)

    project_name = documents.get("project")
    analysis_steps = documents.get("analysis")
    s3_output_files_address = documents.get("upload")
    sample_list = documents.get("sample")

    download_files(project_name, sample_list)

    if "fastqc" in analysis_steps:
        run_fastqc(project_name, sample_list)
    if "bwa-alignment" in analysis_steps:
        run_bwa(project_name, sample_list)

#downloads the data files
def download_files(project_name, sample_list):
    workspace = root_dir + "/scripts/"
    sample_dir = data_dir + "/" + project_name + "/downloads/"
    global data_files_dir
    data_files_dir = sample_dir

    print "downloading files ..."

    for sample_file in sample_list:
        split_files = sample_file.get("filename").split(",")
        for file in split_files:
            file = file.strip()
            subprocess.call(["qsub", workspace + "download.sh", sample_file.get("download"), file, sample_dir])
    PBSTracker.trackPBSQueue(1, "download")

#runs fastqc if in analysis steps
def run_fastqc(project_name, sample_list):
    workspace = root_dir + "/scripts/"
    sample_dir = data_dir + "/" + project_name + "/FastQC/"

    print "executing fastqc..."

    for sample_file in sample_list:
        split_files = sample_file.get("filename").split(",")
        for file in split_files:
            file = file.strip()
            output_file = os.path.splitext(file)[0] + "_fastqc.zip"
            while "." in file:
                file = os.path.splitext(file)[0]
            file += ".fq"
            subprocess.call(["qsub", workspace + "fastqc.sh", data_files_dir + file, sample_dir])
    PBSTracker.trackPBSQueue(1, "fastqc")

#runs bwa if in analysis steps
def run_bwa(project_name, sample_list):
    workspace = root_dir + "/scripts/"
    sample_dir = data_dir + "/" + project_name + "/"

    print "executing bwa alignment ..."

    for sample_file in sample_list:
        split_files = sample_file.get("filename").split(",")
        for file in split_files:
            file = file.strip()
            #subprocess.call(["qsub",    ])

if __name__ == "__main__":
    run_analysis(sys.argv[1])
