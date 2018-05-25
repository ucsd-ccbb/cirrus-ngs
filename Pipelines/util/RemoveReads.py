import sys
import os.path

def parse_inputs(fastq_file):
    root,ext = os.path.splitext(fastq_file)
    return root + ".removed" + ext, root + ".kept" + ext

def remove_reads(fastq_file, removed_reads_file, kept_reads_file):
    f = open(fastq_file, "r")
    g = open(removed_reads_file, "w+")
    h = open(kept_reads_file, "w+")

    prev_header = f.readline()
    line = f.readline()

    while line:
        if len(line.rstrip()) < 18 or len(line.rstrip()) > 25:
            g.write(prev_header)
            g.write(line)
            g.write(f.readline())
            g.write(f.readline())
        else:
            h.write(prev_header)
            h.write(line)
            h.write(f.readline())
            h.write(f.readline())

        prev_header = f.readline()
        line = f.readline()

    f.close()
    g.close()
    h.close()

if __name__ == "__main__":
    fastq_file = sys.argv[1]
    removed_reads_file, kept_reads_file = parse_inputs(fastq_file)
    remove_reads(fastq_file, removed_reads_file, kept_reads_file)
