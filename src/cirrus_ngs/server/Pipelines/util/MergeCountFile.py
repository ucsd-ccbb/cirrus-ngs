__author__ = 'guorongxu'

# This file is used for the merge-counts step in RNA-seq, in the workflow star_htseq and kallisto
import re
import sys
import GencodeGTFParser

def merge_all_sample_count(workflow, sample_list, workspace):

    gencode_gtf_file = "/shared/workspace/software/gencode/gencode.v19.annotation.gtf"
    output_file = workspace + "/all_gene_counts.txt"

    gene_table = GencodeGTFParser.parse(gencode_gtf_file)

    all_gene_counts = []
    samples = get_sample_list(sample_list)  # this is a list of all samples names
    for sample_index, sample_file in enumerate(samples):
        line_index = 0
        #count_file = sample_file.replace(".fastq", "_counts.txt")
        count_file = sample_file + "_counts.txt"

        with open(workspace + "/" + count_file, 'r+') as f:
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

# parse sample string to a list, e.g. "A1 A2 A3" into [A1, A2, A3]
def get_sample_list(sample_list):

    samples = sample_list.split()
    return samples  # return a list of all samples


if __name__ == "__main__":

    workflow = sys.argv[1]
    sample_list = sys.argv[2]
    workspace = sys.argv[3]

    merge_all_sample_count(workflow, sample_list, workspace)
