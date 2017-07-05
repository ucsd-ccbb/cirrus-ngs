__author__ = 'guorongxu'

import re
import sys
import GencodeGTFParser

def merge_all_sample_count(workflow, project_name, sample_list):
    gencode_gtf_file = "/shared/workspace/software/gencode/gencode.v19.annotation.gtf"
    localpath = "/shared/workspace/data_archive/RNASeq/" + project_name + "/" + workflow + "/"
    output_file = localpath + "all_gene_counts.txt"

    gene_table = GencodeGTFParser.parse(gencode_gtf_file)

    all_gene_counts = []
    samples = get_sample_list(sample_list)
    for sample_index, sample_file in enumerate(samples):
        line_index = 0
        count_file = sample_file.replace(".fastq", "_counts.txt")

        with open(localpath + count_file, 'r+') as f:
            lines = f.readlines()
            for line in lines:
                if workflow == "kallisto_deseq_workflow":
                    fields = re.split(r'\t+', line)
                    if fields[0].startswith("gene") or len(fields[0]) == 0:
                        continue

                    if sample_index == 0:
                        all_gene_counts.insert(line_index, [fields[0], fields[1], fields[2], fields[3][:-1]])
                    else:
                        sample_count = all_gene_counts[line_index]
                        sample_count.extend([fields[3][:-1]])
                    line_index = line_index + 1

                elif workflow == "star_htseq_workflow":
                    fields = re.split(r'\t+', line)
                    if sample_index == 0:

                        if fields[0] in gene_table:
                            all_gene_counts.insert(line_index, [gene_table.get(fields[0]), fields[1][:-1]])
                        else:
                            all_gene_counts.insert(line_index, [fields[0], fields[1][:-1]])
                    else:
                        sample_count = all_gene_counts[line_index]
                        sample_count.extend([fields[1][:-1]])
                    line_index = line_index + 1

    filewriter = open(output_file, "a")

    if workflow == "kallisto_deseq_workflow":
        header = "gene\tsymbol\tdescription\t"
    elif workflow == "star_htseq_workflow":
        header = "gene\t"

    for sample_file in samples:
        header = header + sample_file + "\t"
    filewriter.write(header[:-1] + "\n")

    for items in all_gene_counts:
        item_str = ""
        for item in items:
            item_str = item_str + item + "\t"
        filewriter.write(item_str[:-1] + "\n")

    filewriter.close()

def get_sample_list(sample_list):

    samples = []
    fields = re.split(r':+', sample_list)
    for field in fields:
        samples.append(field)

    return samples

if __name__ == "__main__":

    workflow = sys.argv[1]
    project_name = sys.argv[2]
    sample_list = sys.argv[3]

    merge_all_sample_count(workflow, project_name, sample_list)
