__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

def parse(log_file):

    filewriter = open(log_file[:-3] + ".mapping.txt", "w")

    print "system is processing " + log_file
    with open(log_file, 'r+') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("#  novoalign -d"):
                filewriter.write(line[line.find("050316_Wynnis_Tom_mirseq") + len("050316_Wynnis_Tom_mirseq"):line.find("-a TGGAATTCTCGGGTGCCAAGG")] + "\t")
            if line.startswith("#     Read Sequences:  "):
                filewriter.write(line[len("#     Read Sequences:  "):-1] + "\t")
            if line.startswith("#            Aligned:  "):
                filewriter.write(line[len("#            Aligned:  "):-1] + "\n")

    filewriter.close()

if __name__ == "__main__":
    novoalign_log_file = "/Users/guorongxu/Desktop/workspace/NGSProjects/SmallRNASeq/050316_Wynnis_Tom_mirseq/sam/novoalign.log"
    parse(novoalign_log_file)