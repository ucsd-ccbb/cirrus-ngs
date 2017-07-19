__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import sys

def deduplicate(raw_file):
    fastq_id_list = {}

    print "system is processing " + raw_file
    filewriter = open(raw_file.replace("unfilter.fastq", "trim.fastq"), "w")

    with open(raw_file, 'r+') as f:
        lines = f.readlines()
        line_index = 0
        printable = True

        for line in lines:
            if line_index % 4 == 0 and line.rstrip() in fastq_id_list:
                printable = False

            if line_index % 4 == 0 and line.rstrip() not in fastq_id_list:
                fastq_id_list.update({line.rstrip():line.rstrip()})

            if printable:
                filewriter.write(line)

            line_index = line_index + 1

    filewriter.close()
    f.close()

if __name__ == "__main__":

    raw_file = sys.argv[1]
    deduplicate(raw_file)


