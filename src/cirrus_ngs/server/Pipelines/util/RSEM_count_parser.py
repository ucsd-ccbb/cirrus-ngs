__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import re
import os
import sys

def read_data_file(input_file):
    expression_list = {}
    with open(input_file) as f:
        for line_num, line in enumerate(f):
            if line_num == 0:
                expression_list.update({"alignment_statistics":line.rstrip()})
            if line_num == 1:
                expression_list.update({"alignment_certainty":line.rstrip()})
            if line_num == 2:
                expression_list.update({"alignment_Hits":line.rstrip()})
            if line_num > 2:
                # split the string based on tab
                fields = re.split(r'\t+', line)
                expression_list.update({fields[0]:fields[1].rstrip()})

    return expression_list

if __name__ == "__main__":
    workspace = sys.argv[1]
    sample_list = {}
    filewriter = open(workspace + "/all_counts_results.txt", "w")

    for dirpath, directories, filenames in os.walk(workspace):
        for filename in filenames:
            if filename.endswith(".cnt"):       # output of calculate expression
                input_file = os.path.join(dirpath, filename)
                expression_list = read_data_file(input_file)
                sample_list.update({filename:expression_list})
    # no files are found
    if len(sample_list) == 0:
        raise FileNotFoundError("ERROR! No \"cnt\" files are found.")

    filewriter.write("item")
    # sample is filename.cnt
    for sample in sample_list:
        filewriter.write("\t" + sample + "_counts")
    filewriter.write("\n")

    for header in ["alignment_statistics", "alignment_certainty", "alignment_Hits"]:
        filewriter.write(header)
        for sample in sample_list:
            expression_list = sample_list.get(sample)
            expression_values = expression_list.get(header)
            filewriter.write("\t" + expression_values)
        filewriter.write("\n")

    for line_num in range(0, 100):
        filewriter.write(str(line_num) + "\t")
        for sample in sample_list:
            expression_list = sample_list.get(sample)
            if str(line_num) in expression_list:
                filewriter.write(expression_list.get(str(line_num)) + "\t")
            else:
                filewriter.write(str(0) + "\t")
        filewriter.write("\n")

    for header in ["Inf"]:
        filewriter.write(header)
        for sample in sample_list:
            expression_list = sample_list.get(sample)
            expression_values = expression_list.get(header)
            filewriter.write("\t" + expression_values)
        filewriter.write("\n")
