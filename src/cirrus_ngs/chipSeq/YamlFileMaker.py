__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import time

## make a yaml file for analysis of WGSPipeline
def make_yaml_file(yaml_file, workflow, project_name, analysis_steps, s3_input_files_address,
                   sample_list, group_list, s3_output_files_address, genome, style):
    filewriter = open(yaml_file, "w")
    filewriter.write("project: " + project_name + "\n")
    filewriter.write("workflow: " + workflow + "\n")

    analysis_string = "analysis: "
    for analysis in analysis_steps:
        analysis_string = analysis_string + analysis + ", "
    filewriter.write(analysis_string[:-2] + "\n")

    filewriter.write("date: " + time.strftime("%Y/%m/%d") + "\n")
    filewriter.write("upload: " + s3_output_files_address + "\n")
    filewriter.write("genome: " + genome + "\n")
    filewriter.write("style: " + style + "\n")

    for index, group in enumerate(group_list):
        filewriter.write("---\n")
        filewriter.write("sample:\n")

        filename_string = "  filename: "

        if len(group) == 1:
            filewriter.write(filename_string + group[0].rstrip() + "\n")
            filewriter.write("  download: " + s3_input_files_address + "\n")
            filewriter.write("  description: " + group[0].rstrip() + "\n")
            filewriter.write("  group: " + group_list[index] + "\n")
        elif len(group) == 2:
            for sample_name in group:
                filename_string = filename_string + sample_name + ", "
            filewriter.write(filename_string[:-2] + "\n")
            filewriter.write("  download: " + s3_input_files_address + "\n")
            filewriter.write("  description: " + group[0] + "\n")
            filewriter.write("  group: " + str(index) + "\n")


    filewriter.close()


