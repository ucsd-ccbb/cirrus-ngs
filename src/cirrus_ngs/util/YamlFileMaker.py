__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import time

## make a yaml file for analysis of WGSPipeline
def make_yaml_file(yaml_file, project_name, analysis_steps, s3_input_files_address,
                   sample_list, group_list, s3_output_files_address, genome, style):
    filewriter = open(yaml_file, "w")
    filewriter.write("project: " + project_name + "\n")

    analysis_string = "analysis: "
    for analysis in analysis_steps:
        analysis_string = analysis_string + analysis + ", "
    filewriter.write(analysis_string[:-2] + "\n")

    filewriter.write("date: " + time.strftime("%Y/%m/%d") + "\n")
    filewriter.write("upload: " + s3_output_files_address + "\n")
    filewriter.write("genome: " + genome + "\n")
    filewriter.write("style: " + style + "\n")

    for index, sample in enumerate(sample_list):
        filewriter.write("---\n")
        filewriter.write("sample:\n")

        filename_string = "  filename: "

        if len(sample) == 1:
            filewriter.write(filename_string + sample[0].rstrip() + "\n")
            filewriter.write("  download: " + s3_input_files_address + "\n")
            filewriter.write("  description: " + sample[0].rstrip() + "\n")
            filewriter.write("  group: " + group_list[index] + "\n")
        elif len(sample) == 2:
            for sample_name in sample:
                filename_string = filename_string + sample_name + ", "
            filewriter.write(filename_string[:-2] + "\n")
            filewriter.write("  download: " + s3_input_files_address + "\n")
            filewriter.write("  description: " + sample[0].split(".")[0] + "\n")
            filewriter.write("  group: " + group_list[index] + "\n")


    filewriter.close()

def make_wgs_yaml_files(yaml_file, project_name, analysis_steps, 
        s3_input_files_address, sample_list, group_list, 
        s3_output_files_address, genome, style):

    project_string = "project: " + project_name + "\n"
    analysis_string = "analysis: " + ", ".join(analysis_steps) + "\n"
    date_string = "date: " + time.strftime("%Y/%m/%d") + "\n"
    upload_string = "upload: " + s3_output_files_address + "\n"
    genome_string = "genome: " + genome + "\n"
    style_string = "style: " + style + "\n"

    file_names = []

    for index, sample in enumerate(sample_list):
        with open(yaml_file.replace(".yaml", "_{}.yaml".format(group_list[index])), "w") as yaml:
            file_names.append(yaml.name)
            yaml.write(project_string)
            yaml.write(analysis_string)
            yaml.write(date_string)
            yaml.write(upload_string)
            yaml.write(genome_string)
            yaml.write(style_string)
            yaml.write("---\n")
            yaml.write("sample:\n")
            yaml.write("  filename: " + ", ".join(sample) + "\n")
            yaml.write("  download: " + s3_input_files_address + "\n")
            yaml.write("  description: " + sample[0].split(".")[0] + "\n")
            yaml.write("  group: " + group_list[index] + "\n")

    return file_names
