import sys
import os
import subprocess
import PBSTracker
import YamlFileReader

ROOT_DIR = "/shared/workspace/WGSPipeline"
PROJECT_DIR = "/scratch/{}"
DATA_DIR = PROJECT_DIR + "/{}"


##run the actual pipeline
def run_analysis(yaml_file):
    documents = YamlFileReader.parse_yaml_file(yaml_file)

    project_name = documents.get("project")
    analysis_steps = documents.get("analysis")
    s3_output_files_address = documents.get("upload")
    sample_list = documents.get("sample")
    
    global PROJECT_DIR
    PROJECT_DIR = PROJECT_DIR.format(project_name)

    logs_dir = ROOT_DIR.replace("WGSPipeline", "data_archive/DNASeq")
    logs_dir += "/{}/logs".format(project_name)
    
    if not os.path.isdir(logs_dir):
        os.makedirs(logs_dir)

    download_files(project_name, sample_list, logs_dir)

    if "fastqc" in analysis_steps:
        run_fastqc(project_name, sample_list)
    if "bwa-alignment" in analysis_steps:
        run_bwa(project_name, sample_list)

#downloads the data files
def download_files(project_name, sample_list, logs_dir):
    workspace = ROOT_DIR + "/scripts/"

    #print("downloading files ...")

    for sample_file in sample_list:
        split_files = [curr_file.strip() for curr_file in 
                sample_file.get("filename").split(",")]

        global DATA_DIR
        DATA_DIR = DATA_DIR.format(split_files[0].split(".")[0])
        
        if not os.path.isdir(DATA_DIR):
            os.makedirs(DATA_DIR)

        for curr_file in split_files:
            subprocess.call(["qsub", "-N", "download_" + curr_file.split(".")[0],
                workspace + "download.sh", sample_file.get("download"), 
                curr_file, DATA_DIR, logs_dir])

    PBSTracker.trackPBSQueue(1, "download")

#runs fastqc if in analysis steps
def run_fastqc(project_name, sample_list):
    workspace = root_dir + "/scripts/"
    sample_dir = data_dir + "/" + project_name + "/FastQC/"
    global data_files_dir

    print "executing fastqc..."

    for sample_file in sample_list:
        split_files = sample_file.get("filename").split(",")
        for curr_file in split_files:
            curr_file = curr_file.strip()
            output_file = os.path.splitext(curr_file)[0] + "_fastqc.zip"
            while "." in curr_file:
                curr_file = os.path.splitext(curr_file)[0]
            curr_file += ".fq"
            
            global out_log, err_log
            out = open(out_log, "w+")
            err = open(err_log, "w+")
            out.close()
            err.close()

            subprocess.call(["qsub", "-N", curr_file.replace(".", "-") + "_fastqc", 
                "-o", out_log, "-e", err_log, workspace + "fastqc.sh", 
                data_files_dir + curr_file, sample_dir])

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
