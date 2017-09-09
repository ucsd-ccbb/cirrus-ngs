import os
import sys


# abundance file: Path to the abundance file from alignment
def count(sample, abundance_file, output):

    # read_entrez gene
    entrez_gene = os.environ["entrez_gene"]
    genes_by_transcript = {}
    counts_by_gene_id = {}
    # open entrez gene file
    with open(entrez_gene, "r") as f:
        for line in f:
            # remove new line from the end of a line, and split the line by tab
            transcript_id, gene_id = line.rstrip().split("\t")
            genes_by_transcript[transcript_id] = gene_id
            counts_by_gene_id[gene_id] = 0

    # read gene_symbol_description
    gene_description = os.environ["gene_description"]
    gene_2_describe = {}
    # open gene symbol description
    with open(gene_description, "r") as f:
        for line in f:
            gene, symbol, description = line.rstrip().split("\t")
            gene_2_describe[gene] = symbol + "\t" + description

    # read_abundance file
    with open(abundance_file, "r") as f:
        next(f)  # skip header
        for line in f:
            transcript = line.split("\t")[0].split("|")[0]
            counts = line.split("\t")[3]
            counts_by_gene_id[genes_by_transcript[transcript]] += counts

## TODO: translate~~


if __name__ == '__main__':
    count(sys.argv[1])