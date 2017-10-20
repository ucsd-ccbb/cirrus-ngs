import csv
import sys

# This file is for counting reads in miRNA-seq

def count(workspace, sample_name):
    samfile = workspace + "/" + sample_name + ".sam"
    print(samfile)
    sam = open(samfile)

    # initialize the map with all miRNAs
    miRNAmap = {}
    for line in sam.readlines():
        if '@SQ' in line:
            miRNAinfo = line.split('\t')
            # getting the miRNA name (reference sequence name)
            miRNA = miRNAinfo[1].split("SN:", 1)[1]
            miRNAmap[miRNA] = [0]   # initializing the list for that miRNA as 0
            miRNAmap[miRNA].insert(0, miRNA)   # store miRNA name at beginning of the list

    # fill in the mapped reads for each miRNA
    for line in sam.readlines():
        # skip header
        if '@' in line:
            continue
        # get the miRNA mapped/unmapped
        alignmentinfo = line.split('\t')  # split by tab
        miRNA = alignmentinfo[2]

        # 1. if read not mapped
        if miRNA == "*":
            continue

        # if miRNA not in hash map, should never be the case
        if miRNA not in miRNAmap:
            print("ERROR: " + miRNA + " not in samfile!")

        # 2. if read mapped to valid RNA
        else:
            miRNAmap[miRNA][1] += 1  # number of mapped reads to specific miRNA

    sam.close()

    # display counts to output file

    # make a list for the sample
    sample = []
    sample.append(sample_name)

    # create new file and write to it
    with open(workspace + "/" + sample_name + ".counts.out", "w+") as outfile:
        # display header to file
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow('         ',sample)
        for key in miRNAmap:
            # write the value for the key
            writer.writerow(miRNAmap[key])


if __name__ == "__main__":
    count(sys.argv[1], sys.argv[2])
