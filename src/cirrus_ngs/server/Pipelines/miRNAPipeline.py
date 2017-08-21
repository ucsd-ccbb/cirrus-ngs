__author__ = 'Mengyi Liu<mel097@ucsd.edu>'
# this file is in: /shared/workspace/SmallRNASeqPipeline

import sys
import subprocess
import PBSTracker
import YamlFileReader

root = "/scratch"
workspace = "/shared/workspace/SmallRNASeqPipeline/"
scripts = "/shared/workspace/SmallRNASeqPipeline/scripts/"
log = "/shared/workspace/logs/SmallRNASeq"
num_threads = "1"   # append to argument list
min_len = "27"
zipped = "False"


def run_analysis(yaml_file):

    documents = YamlFileReader.parse_yaml_file(yaml_file)

    project_name = documents.get("project")
    analysis_steps = documents.get("analysis")
    s3_input_address = documents.get("upload")
    sample_list = documents.get("sample")

    # makes logs directory
    global log
    log += '/' + project_name

    # call run_fastqc
    if "fastqc" in analysis_steps:
        for sample in sample_list:
            sample_info = get_sample_info(sample)
            print("SAMPLE INFO: ")
            print(sample_info)
            print("project name: " + project_name)
            run_fastqc(project_name, sample_info.get("file_suffix"), root, sample_info.get("fastq_end1"),
                       sample_info.get("fastq_end2"), s3_input_address,
                       sample_info.get("s3_output_address"), log, sample_info.get("is_zipped"))

    # call run_trimmomatic
    if "trim" in analysis_steps:
        for sample in sample_list:
            sample_info = get_sample_info(sample)
            print("SAMPLE INFO: ")
            print(sample_info)
            print("zipped: "+zipped)
            print("project name: " + project_name)
            print("threads: " + num_threads)
            print("min length: "+min_len)
            run_trimmomatic(project_name, sample_info.get("file_suffix"), root, sample_info.get("fastq_end1"),
                            sample_info.get("fastq_end2"), s3_input_address,
                            sample_info.get("s3_output_address"), log, zipped, num_threads, min_len)

    # call cut adapt
    if "cut_adapt" in analysis_steps:
        for sample in sample_list:
            sample_info = get_sample_info(sample)
            print("SAMPLE INFO: ")
            print(sample_info)
            print("zipped: " + zipped)
            print("project name: " + project_name)
            run_cut_adapt(project_name, sample_info.get("file_suffix"), root, sample_info.get("fastq_end1"),
                          sample_info.get("fastq_end2"), s3_input_address,
                          sample_info.get("s3_output_address"), log, zipped)

    # call bowtie2
    if "bowtie2" in analysis_steps:
        for sample in sample_list:
            sample_info = get_sample_info(sample)
            print("SAMPLE INFO: ")
            print(sample_info)
            print("zipped: " + zipped)
            print("project name: " + project_name)
            run_bowtie2(project_name, sample_info.get("file_suffix"), root, sample_info.get("fastq_end1"),
                        sample_info.get("fastq_end2"), s3_input_address,
                        sample_info.get("s3_output_address"), log, zipped)


# process samples
def get_sample_info(sample):
    files = sample.get("filename")
    # set file suffix
    file_suffix = ".fastq"
    if files.find(".fq") > -1:
        file_suffix = ".fq"
    # check if the file is zipped
    is_zipped = "False"
    if files.find(".gz") > -1:
        is_zipped = "True"

    # get forward and reverse end names without suffix
    # fastq_end1 is the forward reads, fastq_end2 is the reverse reads
    file_list = [x.strip() for x in files.split(',')]
    if len(file_list) < 2:
        fastq_end1 = file_list[0][:file_list[0].find(file_suffix)]
        fastq_end2 = "NULL"
    else:
        fastq_end1 = file_list[0][:file_list[0].find(file_suffix)]
        fastq_end2 = file_list[1][:file_list[1].find(file_suffix)]

    # get output address
    s3_output_address = sample.get("download")

    sample_info = {'file_suffix': file_suffix, 'fastq_end1': fastq_end1, 'fastq_end2': fastq_end2,
                   's3_output_address': s3_output_address, 'is_zipped': is_zipped}

    return sample_info


# run fastqc
def run_fastqc(project_name, file_suffix, root_dir, fastq_end1, fastq_end2, s3_input_address, s3_output_address,
               log_dir, is_zipped):

    print("executing fastqc...")
    print(project_name, isinstance(log_dir, str))
    print(file_suffix, isinstance(file_suffix, str))
    print(root_dir, isinstance(root_dir, str))
    print(fastq_end1, isinstance(fastq_end1, str))
    print(fastq_end2, isinstance(fastq_end2, str))
    print(s3_input_address, isinstance(s3_input_address, str))
    print(s3_output_address, isinstance(s3_output_address, str))
    print(log_dir, isinstance(log_dir, str))
    print(is_zipped, isinstance(is_zipped, str))

    subprocess.call(['qsub', "-o", "/dev/null", "-e", "/dev/null", scripts + 'fastqc.sh',
                    project_name, file_suffix, root_dir, fastq_end1, fastq_end2, s3_input_address, s3_output_address,
                    log_dir, is_zipped])

    PBSTracker.trackPBSQueue(1, "fastqc")


# run trimmomatic
def run_trimmomatic(project_name, file_suffix, root_dir, fastq_end1, fastq_end2, s3_input_address, s3_output_address,
                    log_dir, is_zipped, threads, min_length):
    print("executing trimmomatic...")
    print(project_name, isinstance(log_dir, str))
    print(file_suffix, isinstance(file_suffix, str))
    print(root_dir, isinstance(root_dir, str))
    print(fastq_end1, isinstance(fastq_end1, str))
    print(fastq_end2, isinstance(fastq_end2, str))
    print(s3_input_address, isinstance(s3_input_address, str))
    print(s3_output_address, isinstance(s3_output_address, str))
    print(log_dir, isinstance(log_dir, str))
    print(is_zipped, isinstance(is_zipped, str))

    print("threads: "+threads)
    print("min length: "+min_length)

    subprocess.call(['qsub', "-o", "/dev/null", "-e", "/dev/null", scripts + 'trim.sh',
                    project_name, file_suffix, root_dir, fastq_end1, fastq_end2, s3_input_address, s3_output_address,
                    log_dir, is_zipped, threads, min_length])

    PBSTracker.trackPBSQueue(1, "trim")


# run cut adapt
def run_cut_adapt(project_name, file_suffix, root_dir, fastq_end1, fastq_end2, s3_input_address, s3_output_address,
                  log_dir, is_zipped):
    print("executing cut adapt...")
    print(project_name, isinstance(log_dir, str))
    print(file_suffix, isinstance(file_suffix, str))
    print(root_dir, isinstance(root_dir, str))
    print(fastq_end1, isinstance(fastq_end1, str))
    print(fastq_end2, isinstance(fastq_end2, str))
    print(s3_input_address, isinstance(s3_input_address, str))
    print(s3_output_address, isinstance(s3_output_address, str))
    print(log_dir, isinstance(log_dir, str))
    print(is_zipped, isinstance(is_zipped, str))

    subprocess.call(['qsub', "-o", "/dev/null", "-e", "/dev/null", scripts + 'cutadapt.sh',
                    project_name, file_suffix, root_dir, fastq_end1, fastq_end2, s3_input_address, s3_output_address,
                    log_dir, is_zipped])

    PBSTracker.trackPBSQueue(1, "cutadapt")


# run bowtie2 alignment
def run_bowtie2(project_name, file_suffix, root_dir, fastq_end1, fastq_end2, s3_input_address, s3_output_address,
                log_dir, is_zipped):
    print("executing bowtie 2...")
    print(project_name, isinstance(log_dir, str))
    print(file_suffix, isinstance(file_suffix, str))
    print(root_dir, isinstance(root_dir, str))
    print(fastq_end1, isinstance(fastq_end1, str))
    print(fastq_end2, isinstance(fastq_end2, str))
    print(s3_input_address, isinstance(s3_input_address, str))
    print(s3_output_address, isinstance(s3_output_address, str))
    print(log_dir, isinstance(log_dir, str))
    print(is_zipped, isinstance(is_zipped, str))

    subprocess.call(['qsub', "-o", "/dev/null", "-e", "/dev/null", scripts + 'bowtie2.sh',
                    project_name, file_suffix, root_dir, fastq_end1, fastq_end2, s3_input_address, s3_output_address,
                    log_dir, is_zipped])

    PBSTracker.trackPBSQueue(1, "bowtie2")


if __name__ == "__main__":
    run_analysis(sys.argv[1])
