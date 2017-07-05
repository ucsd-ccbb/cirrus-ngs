__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import zipfile
import re
import os

def read_data_file(workspace, file_name):

    archive = zipfile.ZipFile(workspace + "/" + file_name, 'r')
    data_file = archive.read(file_name[:-4] + "/fastqc_data.txt")

    percentages = 0.0

    lines = re.split(r'\n+', data_file)
    printable = False
    for line in lines:
        if line.find("Overrepresented sequences") > -1:
            printable = True
            continue

        if printable and line.find(">>END_MODULE") > -1:
            printable = False

        if printable:
            if line.startswith("#"):
                continue
            else:
                fields = re.split(r'\t+', line)
                if float(fields[2]) > 0.1:
                    percentages = percentages + float(fields[2])

    print file_name + "\t" + str(percentages)

if __name__ == "__main__":
    workspace = "/Users/guorongxu/Desktop/workspace/NGSProjects/SmallRNASeq/042716_Ying_Olefsky_mirseq"

    for root, dirs, files in os.walk(workspace):
        for file in files:
            if file.endswith(".zip"):
                read_data_file(workspace, file)
