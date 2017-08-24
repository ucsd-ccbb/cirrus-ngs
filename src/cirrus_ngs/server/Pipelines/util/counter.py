
# Place this file under /shared/workspace/Pipelines/util

import os
import csv
import sys


# samfiles in "AD2 AD26" format
def count(fafile, countsfilepath, mappingratesfilepath, samfiles, workspace):
    # turn the string into a list, append .sam
    samfiles = list(map(lambda x: x+ ".sam", samfiles.split()))

    numfiles = (len(samfiles))

    fa = open(fafile)

    # intialize map with all miRNAs
    miRNAmap = {}
    for line in fa.readlines():
        # read only name
        if '>' in line:
            RNAinfo = line.split(' ')
            RNA = RNAinfo[0].replace(">", '')
            miRNAmap[RNA] = [0]*(numfiles)
            miRNAmap[RNA].insert(0, RNA)  # store RNA name at beginning of list
            #print RNA

    # init map to count mapping rate for each sample
    mappingRates = []

    # list of sample files
    samples = []

    for i in range(0, numfiles):
        samplename = samfiles[i].replace('.sam', '')
        samples.append(samplename)

        # all samfiles downloaded to workspace
        samfile = open(workspace+"/"+samfiles[i])

        unmapped = 0  # count number of reads not mapped to RNA
        mapped = 0  # count number of reads mapped

        for line in samfile.readlines():
            # skip header
            if '@' in line:
                continue

            # grab RNA
            alignmentinfo = line.split('\t')  # split by tab
            RNA = alignmentinfo[2]

            # 1. if read not mapped
            if RNA == "*":
                unmapped += 1
                continue

            # if miRNA not in hash map, should never be the case
            if RNA not in miRNAmap:
                print ("ERROR: " + RNA + " not in fa file: ")

            # 2. if read mapped to valid RNA
            else:
                mapped += 1  # total number of mapped reads
                miRNAmap[RNA][i+1] += 1  # number of mapped reads to specific RNA

        # add total mapped and unmapped reads to mappingRates map
        # row = sample name, the columns are total number of reads, total number of mapped reads, percentage
        total = mapped + unmapped
        mappingRates.append([samplename, total, mapped, '{0:.1f}%'.format((mapped*100.0)/total)])

        samfile.close()

    # display counts to output file
    with open(countsfilepath, 'w+') as outfile:  # write to output file, + create new file if dne
        # display header to file
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow(samples)
        for key in miRNAmap:
            writer.writerow(miRNAmap[key])

    # display mapping rates to output file
    with open(mappingratesfilepath, 'w+') as outfile:  # write to output file, + create new file if dne
        # display header to file
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow(['     ','Total Reads','Mapped Reads','Mapping Rate'])
        for sample in mappingRates:
            writer.writerow(sample)


if __name__ == "__main__":
    count(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
