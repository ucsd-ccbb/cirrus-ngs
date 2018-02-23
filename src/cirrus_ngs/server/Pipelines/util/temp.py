import sys

#gencode = "/shared/workspace/software/kallisto_util/count_reads/gencode.v23.metadata.EntrezGene"
#gencode = "/shared/workspace/software/references/Hsapiens/hg19/sequence/gencode.v19.metadata.HGNC"
gencode = "/shared/workspace/software/references/Hsapiens/hg19/sequence/gencode.v27lift37.metadata.EntrezGene"
hsa = "/shared/workspace/software/references/Hsapiens/hg19/sequence/new_hsa.txt"
#hsa = "/shared/workspace/software/kallisto_util/count_reads/Hsa_gene_symbol_description.txt"

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
#        sym2descr[sym] = (gene,descr)

sample = sys.argv[1]

#with open("{}.abundance.tsv".format(sample), "r") as f:
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

outfile = "{}_counts1.txt".format(sample)

update_dict = {"100128028":"109729184", "440386":None, "554223":"352962", "83935":"143872", "9142":None}

with open(outfile, "w+") as f:
    f.write("gene\tsymbol\tdescription\t{}\n".format(sample))


    for key in sorted(mash.keys()):
        new_key = update_dict.get(key, key)
        if new_key:
            f.write("{}\t{}\t{}\n".format(key, gene2descr[new_key], mash[key]))
