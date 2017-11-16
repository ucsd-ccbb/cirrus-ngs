import getopt
import sys
from operator import itemgetter

usage_msg = """
usage:
    python SortVCF.py <vcf_file> <chrom_order> [options]

args:
    vcf_file        path to vcf file to be sorted
                        e.g. /scratch/project/raw.vcf
    chrom_order     space-separated string with target chromosome ordering
                        e.g. "1 2 3" or "chr10 chr2 chr5"
                        the target chromsomes must match chromosome names in vcf file

options:
    -o              specify output file, default stdout
    -h              print this usage message and exit
"""

incorrect_chrom_order_error = """
ERROR
    Found chromosome in vcf_file not found in chrom_order input
        chrom_order = \"{}\"
        failing chromosome = \"{}\"
""" + usage_msg


def get_options(in_args):
    optlist, args = getopt.gnu_getopt(in_args, "o:h")
    if "-h" in [x[0] for x in optlist]:
        sys.exit(usage_msg)

    out_file = sys.stdout
    for o, a in optlist:
        if o == "-h":
            sys.exit(usage_msg)
        elif o == "-o":
            out_file = open(a, "w+")

    try:
        vcf_file = args[0]
        chromosome_order = args[1].split()
    except IndexError:
        sys.exit(usage_msg)

    return vcf_file, chromosome_order, out_file

def print_headers(vcf_file, out_fh):
    with open(vcf_file, 'r') as f:
        header_offset = 0
        for line in f:
            if line[0] == "#":
                out_fh.write(line)
                header_offset += len(line)
            else:
                break

    return header_offset

def parse_variants(vcf_file, chromosome_order, header_offset):
    chromosome_dict = {x:[] for x in chromosome_order}
    curr_offset = header_offset
    with open(vcf_file, 'r') as f:
        f.seek(header_offset)
        for line in f:
            entry = line.strip().split()

            try:
                chromosome_dict[entry[0]].append((int(entry[1]), curr_offset))
            except KeyError as err:
                sys.exit(incorrect_chrom_order_error.format(
                    " ".join(chromosome_order), entry[0]))

            curr_offset += len(line)

    return chromosome_dict

def sort_vcf(vcf_file, chromosome_dict, chromosome_order, out_fh):
    with open(vcf_file, "r") as f:
        newline = ""
        for chrom in chromosome_order:
            curr_sorted_variants = sorted(chromosome_dict[chrom], key=itemgetter(0))
            for offset in [x[1] for x in curr_sorted_variants]:
                f.seek(offset)
                out_fh.write(newline + f.readline().rstrip())
                newline = "\n"

if __name__ == "__main__":
    vcf_file, chromosome_order, out_fh = get_options(sys.argv[1:])
    header_offset = print_headers(vcf_file, out_fh)
    c_dict = parse_variants(vcf_file, chromosome_order, header_offset)
    sort_vcf(vcf_file, c_dict, chromosome_order, out_fh)
    out_fh.close()
