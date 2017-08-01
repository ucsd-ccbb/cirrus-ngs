__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import sys
import subprocess
import PBSTracker
import YamlFileReader

root_dir = "/shared/workspace/ChiPSeqPipeline"
data_dir = "/shared/workspace/data_archive/ChiPSeq"

## run all analysis from download, alignment, counting and differential calculation.
def run_analysis(yaml_file):

    documents = YamlFileReader.parse_yaml_file(yaml_file)

    workflow = documents.get("workflow")
    project_name = documents.get("project")
    analysis_steps = documents.get("analysis")
    s3_output_files_address = documents.get("upload")
    style = documents.get("style")
    genome = documents.get("genome")
    sample_list = documents.get("sample")

    ## Download files from s3 and make a design group file.
    download_files(workflow, project_name, sample_list)

    if "fastqc" in analysis_steps:
        run_fastqc(workflow, project_name, sample_list)

    if "alignment" in analysis_steps:
        run_alignment(workflow, project_name, sample_list)

    if "make_tag_directory" in analysis_steps:
        make_tag_directory(workflow, project_name, sample_list)

    if "make_UCSC_file" in analysis_steps:
        make_UCSC_file(workflow, project_name, sample_list)

    if "find_peaks" in analysis_steps:
        find_peaks(workflow, project_name, sample_list, style)

    if "annotate_peaks" in analysis_steps:
        annotate_peaks(workflow, project_name, sample_list, style, genome)

    if "pos2bed" in analysis_steps:
        pos2bed(workflow, project_name, sample_list, style)

    if "find_motifs_genome" in analysis_steps:
        find_motifs_genome(workflow, project_name, sample_list, style, genome)

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

## executing ChipSeq Sequencing alignment
def run_alignment(workflow, project_name, sample_list):
    workspace = root_dir + "/" + workflow + "/scripts/"
    sample_dir = data_dir + "/" + project_name + "/" + workflow + "/"

    print "executing alignment..."
    for sample_file in sample_list:
        if sample_file.get("filename").find(",") > -1:
            sample_1 = sample_file.get("filename")[:sample_file.get("filename").find(",")]
            sample_2 = sample_file.get("filename")[sample_file.get("filename").find(",") + 2:]
            subprocess.call(["qsub", "-pe", "smp", "4", workspace + "alignment.sh", sample_dir + sample_1])
            subprocess.call(["qsub", "-pe", "smp", "4", workspace + "alignment.sh", sample_dir + sample_2])
        else:
            subprocess.call(["qsub", "-pe", "smp", "4", workspace + "alignment.sh", sample_dir + sample_file.get("filename")])
    PBSTracker.trackPBSQueue(1, "alignment")

## make tag directory
def make_tag_directory(workflow, project_name, sample_list):
    workspace = root_dir + "/" + workflow + "/scripts/"
    sample_dir = data_dir + "/" + project_name + "/" + workflow + "/"

    print "executing make_tag_directory..."
    for sample_file in sample_list:
        if sample_file.get("filename").find(",") > -1:
            sample_1 = sample_file.get("filename")[:sample_file.get("filename").find(",")]
            sample_2 = sample_file.get("filename")[sample_file.get("filename").find(",") + 2:]
            input_file = sample_1.replace(".fastq", ".fastq.sam")
            output_tag_folder = sample_1[:-6]
            subprocess.call(["qsub", workspace + "make_tag_directory.sh",
                             sample_dir + output_tag_folder, sample_dir + input_file])
            input_file = sample_2.replace(".fastq", ".fastq.sam")
            output_tag_folder = sample_2[:-6]
            subprocess.call(["qsub", workspace + "make_tag_directory.sh",
                             sample_dir + output_tag_folder, sample_dir + input_file])
        else:
            input_file = sample_file.get("filename").replace(".fastq", ".fastq.sam")
            output_tag_folder = sample_file.get("filename")[:-6]
            subprocess.call(["qsub", workspace + "make_tag_directory.sh",
                             sample_dir + output_tag_folder, sample_dir + input_file])

    PBSTracker.trackPBSQueue(1, "make_tag")

## make UCSC file
def make_UCSC_file(workflow, project_name, sample_list):
    workspace = root_dir + "/" + workflow + "/scripts/"
    sample_dir = data_dir + "/" + project_name + "/" + workflow + "/"

    print "executing make_UCSC_directory..."
    for sample_file in sample_list:
        if sample_file.get("filename").find(",") > -1:
            sample_1 = sample_file.get("filename")[:sample_file.get("filename").find(",")]
            sample_2 = sample_file.get("filename")[sample_file.get("filename").find(",") + 2:]
            input_tag_folder = sample_1[:-6]
            subprocess.call(["qsub", workspace + "make_UCSC_file.sh", sample_dir + input_tag_folder])
            input_tag_folder = sample_2[:-6]
            subprocess.call(["qsub", workspace + "make_UCSC_file.sh", sample_dir + input_tag_folder])
        else:
            input_tag_folder = sample_file.get("filename")[:-6]
            subprocess.call(["qsub", workspace + "make_UCSC_file.sh", sample_dir + input_tag_folder])
    PBSTracker.trackPBSQueue(1, "make_UCSC")

## find peaks
def find_peaks(workflow, project_name, sample_list, style):
    workspace = root_dir + "/" + workflow + "/scripts/"
    sample_dir = data_dir + "/" + project_name + "/" + workflow + "/"

    print "executing find_peaks..."

    for sample_file in sample_list:
        if sample_file.get("filename").find(",") < 0:
            chip_tag_folder = sample_file.get("filename")[:-6]
            if style == "factor":
                output_peak_file = chip_tag_folder + "/peaks.txt"
                subprocess.call(["qsub", workspace + "find_peaks.sh", sample_dir + chip_tag_folder, style,
                        sample_dir + output_peak_file])
            if style == "histone":
                output_peak_file = chip_tag_folder + "/regions.txt"
                subprocess.call(["qsub", workspace + "find_peaks.sh", sample_dir + chip_tag_folder, style,
                        sample_dir + output_peak_file])
        else:
            chip_tag_folder = sample_file.get("filename")[:sample_file.get("filename").find(",")][:-6]
            input_tag_folder = sample_file.get("filename")[sample_file.get("filename").find(",") + 2:][:-6]

            print chip_tag_folder
            print input_tag_folder

            if style == "factor":
                output_peak_file = chip_tag_folder + "_vs_" + input_tag_folder + "/peaks.txt"
                subprocess.call(["qsub", workspace + "find_peaks.sh", sample_dir + chip_tag_folder, style,
                        sample_dir + output_peak_file, sample_dir + input_tag_folder])

            if style == "histone":
                output_peak_file = chip_tag_folder + "_vs_" + input_tag_folder + "/regions.txt"
                subprocess.call(["qsub", workspace + "find_peaks.sh", sample_dir + chip_tag_folder, style,
                        sample_dir + output_peak_file, sample_dir + input_tag_folder])

    PBSTracker.trackPBSQueue(1, "find_peaks")

## annotate peaks
def annotate_peaks(workflow, project_name, sample_list, style, genome):
    workspace = root_dir + "/" + workflow + "/scripts/"
    sample_dir = data_dir + "/" + project_name + "/" + workflow + "/"

    print "executing annotate_peaks..."

    for sample_file in sample_list:
        if sample_file.get("filename").find(",") < 0:
            input_tag_folder = sample_file.get("filename")[:-6]
            if style == "factor":
                subprocess.call(["qsub", workspace + "annotate_peaks.sh", sample_dir + input_tag_folder + "/peaks.txt",
                                 sample_dir + input_tag_folder + "/peaks.annotate.txt", input_tag_folder, genome])

            if style == "histone":
                subprocess.call(["qsub", workspace + "annotate_peaks.sh", sample_dir + input_tag_folder + "/regions.txt",
                                 sample_dir + input_tag_folder + "/regions.annotate.txt", input_tag_folder, genome])
        else:
            chip_tag_folder = sample_file.get("filename")[:sample_file.get("filename").find(",")][:-6]
            input_tag_folder = sample_file.get("filename")[sample_file.get("filename").find(",") + 2:][:-6]
            group_tag_folder = chip_tag_folder + "_vs_" + input_tag_folder
            if style == "factor":
                subprocess.call(["qsub", workspace + "annotate_peaks.sh", sample_dir + group_tag_folder + "/peaks.txt",
                                 sample_dir + group_tag_folder + "/peaks.annotate.txt", sample_dir + group_tag_folder, genome])

            if style == "histone":
                subprocess.call(["qsub", workspace + "annotate_peaks.sh", sample_dir + group_tag_folder + "/regions.txt",
                                 sample_dir + group_tag_folder + "/regions.annotate.txt", sample_dir + group_tag_folder, genome])

    PBSTracker.trackPBSQueue(1, "annotate_p")

## pos2bed
def pos2bed(workflow, project_name, sample_list, style):
    workspace = root_dir + "/" + workflow + "/scripts/"
    sample_dir = data_dir + "/" + project_name + "/" + workflow + "/"

    print "executing pos2bed..."

    for sample_file in sample_list:
        if sample_file.get("filename").find(",") < 0:
            input_tag_folder = sample_file[:-6]
            if style == "factor":
                subprocess.call(["qsub", workspace + "pos2bed.sh", sample_dir + input_tag_folder + "/peaks.txt",
                                 sample_dir + input_tag_folder + "/output.bed"])
            if style == "histone":
                subprocess.call(["qsub", workspace + "pos2bed.sh", sample_dir + input_tag_folder + "/regions.txt",
                                 sample_dir + input_tag_folder + "/output.bed"])
        else:
            chip_tag_folder = sample_file.get("filename")[:sample_file.get("filename").find(",")][:-6]
            input_tag_folder = sample_file.get("filename")[sample_file.get("filename").find(",") + 2:][:-6]
            group_tag_folder = chip_tag_folder + "_vs_" + input_tag_folder
            if style == "factor":
                subprocess.call(["qsub", workspace + "pos2bed.sh", sample_dir + group_tag_folder + "/peaks.txt",
                                 sample_dir + group_tag_folder + "/output.bed"])

            if style == "histone":
                subprocess.call(["qsub", workspace + "pos2bed.sh", sample_dir + group_tag_folder + "/regions.txt",
                                 sample_dir + group_tag_folder + "/output.bed"])

        PBSTracker.trackPBSQueue(1, "pos2bed")

## find motifs genome
def find_motifs_genome(workflow, project_name, sample_list, style, genome):
    workspace = root_dir + "/" + workflow + "/scripts/"
    sample_dir = data_dir + "/" + project_name + "/" + workflow + "/"

    print "executing find_motifs_genome..."

    for sample_file in sample_list:
        if sample_file.get("filename").find(",") < 0:
            input_tag_folder = sample_file.get("filename")[:-6]
            if style == "factor":
                subprocess.call(["qsub", workspace + "find_motifs_genome.sh", sample_dir + input_tag_folder + "/peaks.txt",
                                 sample_dir + input_tag_folder + "/MotifOutput/", genome])

            if style == "histone":
                subprocess.call(["qsub", workspace + "find_motifs_genome.sh", sample_dir + input_tag_folder + "/regions.txt",
                                 sample_dir + input_tag_folder + "/MotifOutput/", genome])
        else:
            chip_tag_folder = sample_file.get("filename")[:sample_file.get("filename").find(",")][:-6]
            input_tag_folder = sample_file.get("filename")[sample_file.get("filename").find(",") + 2:][:-6]
            group_tag_folder = chip_tag_folder + "_vs_" + input_tag_folder
            if style == "factor":
                subprocess.call(["qsub", workspace + "find_motifs_genome.sh", sample_dir + group_tag_folder + "/peaks.txt",
                                 sample_dir + group_tag_folder + "/MotifOutput/", genome])

            if style == "histone":
                subprocess.call(["qsub", workspace + "find_motifs_genome.sh", sample_dir + group_tag_folder + "/regions.txt",
                                 sample_dir + group_tag_folder + "/MotifOutput/", genome])

        PBSTracker.trackPBSQueue(1, "find_motif")

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
