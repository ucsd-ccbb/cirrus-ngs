__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import os

if __name__ == "__main__":

    sam_dir = "/Users/guorongxu/Desktop/workspace/NGSProjects/SmallRNASeq/050316_Wynnis_Tom_mirseq/sam"

    for root, dirs, files in os.walk(sam_dir):
        for file in files:
            if file.endswith(".sam"):
                sam_file = (os.path.join(root, file))
                os.popen("qsub " + sam_dir + "/count.sh " + sam_file)