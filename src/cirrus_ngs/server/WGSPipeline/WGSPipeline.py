import sys
import os
import subprocess
import PBSTracker
import YamlFileReader

ROOT_DIR = "/shared/workspace/WGSPipeline"
PROJECT_DIR = "/scratch/{}"
DATA_DIR = ""


##run the actual pipeline
def run_analysis(yaml_file):
    documents = YamlFileReader.parse_yaml_file(yaml_file)

    project_name = documents.get("project")
    analysis_steps = documents.get("analysis")
    s3_output_files_address = documents.get("upload")
    sample_list = documents.get("sample")
    
    global PROJECT_DIR
    PROJECT_DIR = PROJECT_DIR.format(project_name)

    print(PROJECT_DIR)

    logs_dir = ROOT_DIR.replace("WGSPipeline", "data_archive/DNASeq")
    logs_dir += "/{}/logs".format(project_name)
    
    if not os.path.isdir(logs_dir):
        os.makedirs(logs_dir)

    download_files(project_name, sample_list, logs_dir)

    if "fastqc" in analysis_steps:
        run_fastqc(project_name, sample_list, logs_dir)
    if not "notrim" in analysis_steps:
        trim(project_name, sample_list, logs_dir)
    if "bwa-alignment" in analysis_steps:
        run_bwa(project_name, sample_list)

#downloads the data files
def download_files(project_name, sample_list, logs_dir):
    workspace = ROOT_DIR + "/scripts/"

    #print("downloading files ...")
    sample_file = sample_list[0]

    split_files = sample_file.get("filename").split(",")

    files = [split_files[x].strip() if x < len(split_files)
            else "NULL" for x in range(2)]
    
    global DATA_DIR, PROJECT_DIR
    DATA_DIR = PROJECT_DIR + "/" + split_files[0].split(".")[0]
    
    if not os.path.isdir(DATA_DIR):
        os.makedirs(DATA_DIR)

    subprocess.call(["bash", workspace + "download.sh", sample_file.get("download"), 
        PROJECT_DIR, files[0], files[1], logs_dir])

    print("Finished downloading {} files".format(sample_list[0].get("group")))

#runs fastqc if in analysis steps
def run_fastqc(project_name, sample_list, logs_dir):
    workspace = ROOT_DIR + "/scripts/"

    sample_file = sample_list[0]

    split_files = sample_file.get("filename").split(",")

    files = [split_files[x].strip() if x < len(split_files)
            else "NULL" for x in range(2)]

    subprocess.call(["bash", workspace + "fastqc.sh", DATA_DIR, PROJECT_DIR, files[0],
        files[1], logs_dir])

    print("Finished performing {} FastQC analysis".format(sample_file.get("group")))

def trim(project_name, sample_list, logs_dir):
    workspace = ROOT_DIR + "/scripts/"

    sample_file = sample_list[0]
    split_files = sample_file.get("filename").split(",")
    files = [split_files[x].strip() if x < len(split_files)
            else "NULL" for x in range(2)]

    subprocess.call(["bash", workspace + "trim.sh", DATA_DIR, PROJECT_DIR,
        files[0], files[1], logs_dir])

    print("Finished trimming {}".format(sample_file.get("group")))

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
