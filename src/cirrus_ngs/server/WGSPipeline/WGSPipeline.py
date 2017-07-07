import sys
import os
import subprocess
import PBSTracker
import YamlFileReader

root_dir = "/shared/workspace/WGSPipeline"
data_dir = "/shared/workspace/data_archive/DNASeq"
data_files_dir = ""
out_log = data_dir + "/%s/logs/$JOB_NAME_$JOB_ID.o"
err_log = data_dir + "/%s/logs/$JOB_NAME_$JOB_ID.e"

##run the actual pipeline
def run_analysis(yaml_file):
    documents = YamlFileReader.parse_yaml_file(yaml_file)

    project_name = documents.get("project")
    analysis_steps = documents.get("analysis")
    s3_output_files_address = documents.get("upload")
    sample_list = documents.get("sample")
    
    logs_dir = data_dir + "/%s/logs/" % project_name
    global out_log
    global err_log
    out_log = out_log % project_name
    err_log = err_log % project_name
    
    if not os.path.isdir(logs_dir):
        os.makedirs(logs_dir)
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

    #print "downloading files ..."

    for sample_file in sample_list:
        split_files = sample_file.get("filename").split(",")
        for curr_file in split_files:
            curr_file = curr_file.strip()

            global out_log, err_log
            out= open(out_log, "w+")
            err= open(err_log, "w+")
            out.close()
            err.close()

            subprocess.call(["qsub", "-N", curr_file.replace(".", "-") + "_download",
                "-o", out_log, "-e", err_log, workspace + "download.sh",
                sample_file.get("download"), curr_file, sample_dir])

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
