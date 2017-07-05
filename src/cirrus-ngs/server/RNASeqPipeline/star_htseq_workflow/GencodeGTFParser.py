__author__ = 'guorongxu'

import re

def parse(gencode_gtf_file):

    print "system is processing " + gencode_gtf_file

    gene_table = {}
    with open(gencode_gtf_file, 'r+') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("##"):
                continue
            fields = re.split(r'\t+', line)

            if fields[2] == "gene":
                gene_id = line[line.find("gene_id") + len("gene_id") + 2:line.find("transcript_id")-3]
                gene_name = line[line.find("gene_name") + len("gene_name") + 2:line.find("transcript_type")-3]

                if gene_id not in gene_table:
                    gene_table.update({gene_id:gene_name})

    return gene_table
