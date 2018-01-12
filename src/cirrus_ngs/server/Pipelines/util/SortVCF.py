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
    """Parses the command line arguments.

    Makes sure that the required arguments
    are present and returns interpretations
    of the required arguments.

    Args:
        in_args: argv list without the 0th element
        
    Returns:
        tuple of form
        (
            string path to vcf_file to be sorted,
            list of string chromosomes that represent the sorted order,
            file handle to output file
        )

    Raises:
        SystemExit if the help flag is used
        SystemExit if there aren't enough command line arguments
    """
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
    """Prints the vcf header to the target output file

    VCF file headers don't need to be sorted, so they
    are written as is to the intended output stream.

    Args:
        vcf_file: string path to vcf file to be sorted
        out_fh: file handle for output file

    Returns: 
        the offset in vcf_file file that corresponds
        to the start of non-header data
    """
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
    """Creates a dictionary that contains file offsets for each chromsome.

    In order to sort by some arbitrary order the VCF file is 
    run through once to store the line offsets for each
    entry for a chromosome in a list. This data structure
    can be sorted without parsing the whole file again, since
    the offsets will lead to the proper information. Actual lines
    aren't stored to prevent excessive memory usage.

    Args:
        vcf_file: string path to vcf_file to be sorted
        chromosome_order: list of string chromsomes in the intended sorting order
        header_offset: the character in the vcf_file that represents the start of non-header data

    Returns:
        dictionary in the form of 
        {
            chromosome: [
                            (int location of current variant in chromsome,
                            int offset in vcf_file of current variant)
                        ]
        }

    Raises:
        SystemExit if chromosome not in chromosome_order is found in vcf_file
    """
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
    """Writes sorted VCF data to the output file.

    Sorting is based on output dictionary of parse_variants function.
    For each chromsome the variants are ordered by their position in 
    that chromosome (first item in tuple in value for chromsome_dict).
    The results are written to the output file in order.

    Args:
        vcf_file: string path to the VCF file to be sorted
        chromsome_dict: dictionary with offset/chromsome data. See parse_variants for structure.
        out_fh: file handle of output file

    Returns:
        None
    """
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
