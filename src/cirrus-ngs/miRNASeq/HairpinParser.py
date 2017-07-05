__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import re

def parse(hairpin_file, output_file):

    print "system is processing " + hairpin_file
    filewriter = open(output_file, "w")

    with open(hairpin_file, 'r+') as f:
        printable = False
        lines = f.readlines()
        for line in lines:
            if line.find(">") > -1:
                printable = False

            if line.find(">") > -1 and line.find("musculus") > -1:
                printable = True
            if printable:
                filewriter.write(line)

    filewriter.close()

if __name__ == "__main__":

    hairpin_file = "/Users/guorongxu/Desktop/jupyter-genomics/data/miRNASeq/hairpin.fa"
    output_file = "/Users/guorongxu/Desktop/jupyter-genomics/data/miRNASeq/hairpin_musculus.fa"
    parse(hairpin_file, output_file)