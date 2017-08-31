__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import time

## make a yaml file for analysis of any Pipeline
def make_yaml_file(yaml_file, project_name, analysis_steps, s3_input_files_address,
        sample_list, group_list, s3_output_files_address, genome, style, pairs_list):
    filewriter = open(yaml_file, "w")
    filewriter.write("project: " + project_name + "\n")

    filewriter.write("analysis:\n")
    for analysis in analysis_steps:
        filewriter.write("  - {}\n".format(analysis))

    filewriter.write("date: " + time.strftime("%Y/%m/%d") + "\n")
    filewriter.write("upload: " + s3_output_files_address + "\n")
    filewriter.write("download: " + s3_input_files_address + "\n")
    filewriter.write("genome: " + genome + "\n")
    filewriter.write("style: " + style + "\n")
    filewriter.write("pairs:\n")
    if not pairs_list:
        filewriter.write("  {}\n")
    else:
        for normal, tumor in pairs_list.items():
            filewriter.write("  {}: {}\n".format(normal, tumor))

    
    filewriter.write("sample:\n")
    for index, sample in enumerate(sample_list):
        filename_string = "- filename: "

        if len(sample) == 1:
            filewriter.write(filename_string + sample[0].rstrip() + "\n")
            filewriter.write("  description: " + sample[0].split(".")[0] + "\n")
            filewriter.write("  group: " + group_list[index] + "\n")
        elif len(sample) == 2:
            for sample_name in sample:
                filename_string = filename_string + sample_name + ", "
            filewriter.write(filename_string[:-2] + "\n")
            filewriter.write("  description: " + sample[0].split(".")[0] + "\n")
            filewriter.write("  group: " + group_list[index] + "\n")

    filewriter.close()
