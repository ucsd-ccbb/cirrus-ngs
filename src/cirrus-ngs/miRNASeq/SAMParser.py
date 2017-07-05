__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import os
import re
import sys

def parse(sam_file):
    mirna_list = {}
    filewriter = open(sam_file[:-4] + ".counts.txt", "w")

    print "system is processing " + sam_file
    with open(sam_file, 'r+') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("@"):
                continue
            else:
                fields = re.split(r'\t+', line)
                if fields[2] != "*" and fields[2] not in mirna_list:
                    mirna_list.update({fields[2]:1})
                elif fields[2] != "*" and fields[2] in mirna_list:
                    count = mirna_list.get(fields[2])
                    mirna_list.update({fields[2]:count + 1})

    filewriter.write("mirna_name\tcount\n")
    for mirna in mirna_list:
        count = mirna_list.get(mirna)
        filewriter.write(mirna + "\t" + str(count) + "\n")

    filewriter.close()

if __name__ == "__main__":
    sam_file = sys.argv[1]
    parse(sam_file)