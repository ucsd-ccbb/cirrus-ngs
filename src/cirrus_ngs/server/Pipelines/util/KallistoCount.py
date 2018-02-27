import sys

gencode = "/shared/workspace/software/references/Hsapiens/hg19/annotations/gencode.v27lift37.metadata.EntrezGene"

# this is just a trimmed version of Homo_sapiens.gene_info from Entrez Gene
hsa = "/shared/workspace/software/references/Hsapiens/hg19/annotations/new_hsa.txt"

genehash = {}
mash = {}

with open(gencode, "r") as f:
    for line in f:
        line = line.rstrip()
        t, g = line.split("\t")[:2]
        genehash[t] = g
        mash[g] = 0

gene2descr = {}
with open(hsa, "r") as f:
    for line in f:
        line = line.rstrip()
        gene, sym, descr = line.split("\t")[:3]
        gene2descr[gene] = "{}\t{}".format(sym, descr)

sample = sys.argv[1]

with open(sample, "r") as f:
    f.readline() # skip header
    for line in f:
        line = line.rstrip()
        line = line.split("\t")
        transcr = line[0]
        count = line[3]
        transcr = transcr.split("|")[0]
        if transcr in genehash:
            mash[genehash[transcr]] += float(count)

outfile = "{}_counts.txt".format(sample.split(".")[0])

# gencode is a little behind entrez gene
# these are the genes that were deprecated/updated 
update_dict = {"100128028":"109729184", "440386":None, "554223":"352962", "83935":"143872", "9142":None}

with open(outfile, "w+") as f:
    f.write("gene\tsymbol\tdescription\t{}\n".format(sample))

    for key in sorted(mash.keys(), key=int):
        new_key = update_dict.get(key, key)
        if new_key:
            f.write("{}\t{}\t{}\n".format(key, gene2descr[new_key], mash[key]))
