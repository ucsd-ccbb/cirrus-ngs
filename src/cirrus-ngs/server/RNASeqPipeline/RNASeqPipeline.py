__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import sys
import subprocess
import PBSTracker
import GroupFileMaker
import YamlFileReader

root_dir = "/shared/workspace/RNASeqPipeline"
data_dir = "/shared/workspace/data_archive/RNASeq"

## run all analysis from download, alignment, counting and differential calculation.
def run_analysis(yaml_file):
    documents = YamlFileReader.parse_yaml_file(yaml_file)
    workflow = documents.get("workflow")
    project_name = documents.get("project")
    analysis_steps = documents.get("analysis")
    s3_output_files_address = documents.get("upload")
    sample_list = documents.get("sample")

    ## Download files from s3 and make a design group file.
    download_files(workflow, project_name, sample_list)
    make_group_file(workflow, project_name, sample_list)

    if "fastqc" in analysis_steps:
        run_fastqc(workflow, project_name, sample_list)

    if "alignment" in analysis_steps:
        ## Before alignment, we need to trim the fastq file first
        fastq_trim(workflow, project_name, sample_list)
        run_alignment(workflow, project_name, sample_list)

    if "variant_calling" in analysis_steps:
        run_variant_calling(workflow, project_name, sample_list)

    if "counting" in analysis_steps:
        count_reads(workflow, project_name, sample_list)
        merge_count_files(workflow, project_name, sample_list)

    if "differential_calculation" in analysis_steps:
        differential_calculate(workflow, project_name)

    ## Upload resulting files to s3.
    upload_files(workflow, project_name, s3_output_files_address)

    print "======================================================"
    print "The processing of the project \"" + project_name + "\" is done!"
    print "======================================================"

## download file from s3
def download_files(workflow, project_name, sample_list):
    workspace = root_dir + "/" + workflow + "/scripts/"
    sample_dir = data_dir + "/" + project_name + "/" + workflow + "/"

    print "executing download files..."

    ##copying data from s3 to local drive
    for sample_file in sample_list:

        if sample_file.get("filename").find(",") > -1:
            sample_1 = sample_file.get("filename")[:sample_file.get("filename").find(",")]
            sample_2 = sample_file.get("filename")[sample_file.get("filename").find(",") + 2:]
            subprocess.call(["qsub", workspace + "download.sh", sample_file.get("download"),
                         sample_1, sample_dir])
            subprocess.call(["qsub", workspace + "download.sh", sample_file.get("download"),
                         sample_2, sample_dir])
        else:
            subprocess.call(["qsub", workspace + "download.sh", sample_file.get("download"),
                         sample_file.get("filename"), sample_dir])
    PBSTracker.trackPBSQueue(1, "download")

## making design file for all samples.
def make_group_file(workflow, project_name, sample_list):
    sample_dir = data_dir + "/" + project_name + "/" + workflow + "/"
    print "making group file..."
    GroupFileMaker.make_group_file(workflow, project_name, sample_dir + "group.txt", sample_list)

    ## viewing the group file.
    subprocess.call(["head", sample_dir + "group.txt"])

## running fastqc for all samples
def run_fastqc(workflow, project_name, sample_list):
    workspace = root_dir + "/" + workflow + "/scripts/"
    sample_dir = data_dir + "/" + project_name + "/" + workflow + "/"

    print "executing fastqc..."
    for sample_file in sample_list:
        if sample_file.get("filename").find(",") > -1:
            sample_1 = sample_file.get("filename")[:sample_file.get("filename").find(",")]
            sample_2 = sample_file.get("filename")[sample_file.get("filename").find(",") + 2:]
            output_file = sample_1.replace(".fastq", "_fastqc.zip")
            subprocess.call(["qsub", workspace + "fastqc.sh", sample_dir + sample_1,
                             sample_dir + output_file])
            output_file = sample_2.replace(".fastq", "_fastqc.zip")
            subprocess.call(["qsub", workspace + "fastqc.sh", sample_dir + sample_2,
                             sample_dir + output_file])
        else:
            output_file = sample_file.get("filename").replace(".fastq", "_fastqc.zip")
            subprocess.call(["qsub", workspace + "fastqc.sh", sample_dir + sample_file.get("filename"),
                             sample_dir + output_file])
    PBSTracker.trackPBSQueue(1, "fastqc")

## executing fastq trimming
def fastq_trim(workflow, project_name, sample_list):
    workspace = root_dir + "/" + workflow + "/scripts/"
    sample_dir = data_dir + "/" + project_name + "/" + workflow + "/"

    print "executing fastq trimming..."
    for sample_file in sample_list:
        if sample_file.get("filename").find(",") > -1:
            fastq_name = sample_file.get("filename")[:sample_file.get("filename").find(",")]
            input_file = fastq_name[:fastq_name.rfind("_")]
            output_file = fastq_name[:fastq_name.rfind("_")]

            subprocess.call(["qsub", workspace + "fastq_trim_pe.sh", sample_dir + input_file,
                             sample_dir + output_file])
        else:
            input_file = sample_file.get("filename").replace(".fastq", "")
            output_file = sample_file.get("filename").replace(".fastq", "")

            subprocess.call(["qsub", workspace + "fastq_trim_se.sh", sample_dir + input_file,
                             sample_dir + output_file])

    PBSTracker.trackPBSQueue(1, "fastq_trim")

## executing RNA Sequencing alignment
def run_alignment(workflow, project_name, sample_list):
    workspace = root_dir + "/" + workflow + "/scripts/"
    sample_dir = data_dir + "/" + project_name + "/" + workflow + "/"

    print "executing alignment..."
    for sample_file in sample_list:
        print "alignment sample file: " + sample_file.get("filename")
        if sample_file.get("filename").find(",") > -1:
            fastq_file = sample_file.get("filename")[:sample_file.get("filename").find(",")]
            input_file = fastq_file[:fastq_file.rfind("_")]
            output_file = fastq_file[:fastq_file.rfind("_")]

            subprocess.call(["qsub", "-pe", "smp", "4", workspace + "alignment_pe.sh", sample_dir + input_file,
                             sample_dir + output_file])

        else:
            input_file = sample_file.get("filename").replace(".fastq", "")
            output_file = sample_file.get("filename").replace(".fastq", "")

            subprocess.call(["qsub", "-pe", "smp", "4", workspace + "alignment_se.sh", sample_dir + input_file,
                             sample_dir + output_file])
    PBSTracker.trackPBSQueue(1, "alignment")

## executing gatk variants calling
def run_variant_calling(workflow, project_name, sample_list):
    workspace = root_dir + "/" + workflow + "/"

    print "executing variant calling..."
    for sample_file in sample_list:
        if sample_file.get("filename").find(",") > -1:
            fastq_file = sample_file.get("filename")[:sample_file.get("filename").find(",")]
            output_file = fastq_file[:fastq_file.rfind("_")]
            subprocess.call(["qsub", workspace + "scripts/variant_calling.sh", output_file, project_name, workflow])
        else:
            output_file = sample_file.get("filename").replace(".fastq", "")
            subprocess.call(["qsub", workspace + "scripts/variant_calling.sh", output_file, project_name, workflow])
    PBSTracker.trackPBSQueue(1, "variant")

## counting reads for each gene or transcript
def count_reads(workflow, project_name, sample_list):
    workspace = root_dir + "/" + workflow + "/"
    sample_dir = data_dir + "/" + project_name + "/" + workflow + "/"

    print "executing counting..."
    for sample_file in sample_list:
        fastq_file = sample_file.get("filename")
        if sample_file.get("filename").find(",") > -1:
            fastq_file = sample_file.get("filename")[:sample_file.get("filename").find(",")]

        if workflow == "star_htseq_workflow":
            input_file = fastq_file.replace(".fastq", "/Aligned.out.sorted.bam")
            output_file = fastq_file.replace(".fastq", "_counts.txt")
            subprocess.call(["qsub", workspace + "scripts/count_reads.sh",
                        sample_dir + input_file, sample_dir + output_file])
        elif workflow == "kallisto_deseq_workflow":
            input_file = fastq_file.replace(".fastq", "")
            output_file = fastq_file.replace(".fastq", "")
            subprocess.call(["qsub", workspace + "scripts/count_reads.sh",
                        sample_dir + input_file, sample_dir + output_file])

    PBSTracker.trackPBSQueue(1, "count_read")

## merging all count files for all samples
def merge_count_files(workflow, project_name, sample_list):
    workspace = root_dir + "/" + workflow + "/"
    sample_dir = data_dir + "/" + project_name + "/" + workflow + "/"

    print "executing merge count files..."

    sample_str = ""
    for sample in sample_list:
        fastq_file = sample.get("filename")
        if sample.get("filename").find(",") > -1:
            fastq_file = sample.get("filename")[:sample.get("filename").find(",")]
        sample_str = sample_str + fastq_file + ":"

    subprocess.call(["qsub", workspace + "scripts/merge_counts.sh", workflow, project_name, sample_str[:-1]])

    PBSTracker.trackPBSQueue(1, "merge_coun")

    ## viewing count files
    subprocess.call(["head", sample_dir + "all_gene_counts.txt"])


## running differential calculation.
def differential_calculate(workflow, project_name):
    workspace = root_dir + "/" + workflow + "/"
    sample_dir = data_dir + "/" + project_name + "/" + workflow + "/"

    print "executing differential calculating..."

    subprocess.call(["qsub", workspace + "/scripts/RNA-seq_limma.sh", sample_dir + "group.txt", sample_dir])

## uploading resulting files to s3.
def upload_files(workflow, project_name, s3_output_files_address):

    print "executing upload files..."

    workspace = root_dir + "/" + workflow + "/scripts/"
    sample_dir = data_dir + "/" + project_name + "/" + workflow + "/"
    subprocess.call(["qsub", workspace + "upload.sh", sample_dir, s3_output_files_address + "/" + project_name + "/" + workflow])

    PBSTracker.trackPBSQueue(1, "upload")

if __name__ == '__main__':
    #yaml_file = "/Users/guorongxu/Desktop/workspace/projects/jupyter-genomics_bitbucket/src/awsCluster/chipSeq/Sample_cDNA.yaml"
    yaml_file = sys.argv[1]
    run_analysis(yaml_file)
