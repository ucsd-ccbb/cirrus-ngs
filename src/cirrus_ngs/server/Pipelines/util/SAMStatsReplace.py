__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import os
import sys

if __name__ == '__main__':
    counts_file = sys.argv[1]
    sam_stats_file = sys.argv[2]
    new_sam_stats_file = sam_stats_file.replace(".txt", ".rep.txt")

    filewriter = open(new_sam_stats_file, "w")

    fastq_counts = ""
    with open(counts_file, 'r+') as f:
        lines = f.readlines()

        for line in lines:
            fastq_counts = line.rstrip()
            break
    f.close()

    with open(sam_stats_file, 'r+') as f:
        lines = f.readlines()

        for line in lines:
            if line.find("raw total sequences") > -1:
                filewriter.write("SN\traw total sequences:\t" + fastq_counts + "\n")
            else:
                filewriter.write(line)

    filewriter.close()

    os.remove(sam_stats_file)
    os.rename(new_sam_stats_file, sam_stats_file)
