__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import os
import re
import sys

def parse(count_file_dir):
    mirna_count_list = {}
    count_file_list = []
    file_name_list = []

    for root, dirs, files in os.walk(count_file_dir):
        for file in files:
            if file.endswith(".counts.txt"):
                count_file_list.append(os.path.join(root, file))
                file_name_list.append(file)

    total_file_num = len(count_file_list)

    for count_index, count_file in enumerate(count_file_list):
        print "system is processing " + count_file
        with open(count_file, 'r+') as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith("mirna_name"):
                    continue
                else:
                    fields = re.split(r'\t+', line)
                    if fields[0] not in mirna_count_list:
                        count_list = []
                        for index in range(0, total_file_num):
                            count_list.append(0)
                        count_list[count_index] = fields[1][:-1]
                        mirna_count_list.update({fields[0]:count_list})
                    elif fields[0] in mirna_count_list:
                        count_list = mirna_count_list.get(fields[0])
                        count_list[count_index] = fields[1][:-1]

    filewriter = open(count_file_dir + "/" + "miRNA.all.counts.txt", "w")
    head_string = "mirna_name\t"
    for file_name in file_name_list:
        head_string = head_string + file_name + "\t"
    filewriter.write(head_string[:-1] + "\n")

    for mirna in mirna_count_list:
        count_string = mirna + "\t"
        count_list = mirna_count_list.get(mirna)
        for count in count_list:
            count_string = count_string + str(count) + "\t"
        filewriter.write(count_string[:-1] + "\n")

    filewriter.close()

if __name__ == "__main__":

    count_file_dir = sys.argv[1]
    parse(count_file_dir)

