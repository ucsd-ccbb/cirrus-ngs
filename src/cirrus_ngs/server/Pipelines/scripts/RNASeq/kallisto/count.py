import sys

entrez_gene = "/shared/workspace/software/kallisto/count_reads/gencode.v23.metadata.EntrezGene"
hsa_gene_symbol_des = "/shared/workspace/software/kallisto/count_reads/Hsa_gene_symbol_description.txt"

def count(input_address):
    read_entrez()
    read_gene_symbol_des()

    # with open(input_address):

def read_entrez():
    genes_by_transcript = {}
    counts_by_gene_id = {}
    # open entrez gene file
    with open(entrez_gene, "r") as f:
        for line in f:
            # TODO remove new line from the end of a line

            transcript_id, gene_id = line.split("\t")
            genes_by_transcript[transcript_id] = gene_id
            counts_by_gene_id[gene_id] = 0


def read_gene_symbol_des():
    gene_2_describe = {}
    # open gene symbol description
    with open(hsa_gene_symbol_des, "r") as f:
        for line in f:
            gene, symbol, description = line.split("\t")
            gene_2_describe[gene] = symbol + "\t" + description


if __name__ == '__main__':
    count(sys.argv[1])