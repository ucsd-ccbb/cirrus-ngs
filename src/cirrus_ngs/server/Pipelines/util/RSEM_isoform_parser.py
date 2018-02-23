__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import re
import os
import sys

def read_data_file(input_file):
    expression_list = []
    with open(input_file) as f:
        for line in f:
            if not line.startswith("transcript_id"):
               expression_list.append(line)
    return expression_list

if __name__ == "__main__":
    workspace = sys.argv[1]
    sample_list = {}
    filewriter = open(workspace + "/all_isoforms_results.txt", "w")

    for dirpath, directories, filenames in os.walk(workspace):
        for filename in filenames:
            if filename.endswith(".isoforms.results"):      # output of calculate expression
                input_file = os.path.join(dirpath, filename)
                expression_list = read_data_file(input_file)
                sample_list.update({filename:expression_list})
    # no files are found
    if len(sample_list) == 0:
        raise FileNotFoundError("ERROR! No \"isoforms.results\" files are found.")

    filewriter.write("transcript_id\tgene_id")
    for sample in sample_list:
        filewriter.write("\t" + sample + "_length\t" + sample + "_effective_length\t" + sample + "_expected_count\t" + sample + "_TPM\t" + sample + "_FPKM\t" + sample + "_IsoPct")
    filewriter.write("\n")

    for line_num in range(0, len(expression_list)):
        for index, sample in enumerate(sample_list):
            expression_list = sample_list.get(sample)
            expression_values = expression_list[line_num]
            fields = re.split(r'\t+', expression_values)
            if index == 0:
                filewriter.write(fields[0] + "\t" + fields[1] + "\t" + fields[2] + "\t" + fields[3] + "\t" + fields[4] + "\t" + fields[5] + "\t" + fields[6] + "\t" + fields[7].rstrip())
            else:
                filewriter.write("\t" + fields[2] + "\t" + fields[3] + "\t" + fields[4] + "\t" + fields[5] + "\t" + fields[6] + "\t" + fields[7].rstrip())

        filewriter.write("\n")
