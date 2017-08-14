__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

def configure_star(option_list):
    output_file_name = "/Users/guorongxu/Desktop/JupyterPythonNotebook/alignment.sh"

    filewriter = open(output_file_name, "w")
    filewriter.write("#!/bin/bash\n")
    filewriter.write("\n")
    filewriter.write("samtools=/shared/workspace/software/samtools/samtools-1.1/samtools\n")
    filewriter.write("STAR=/shared/workspace/software/STAR/2.5.1a/STAR/bin/Linux_x86_64/STAR\n")
    filewriter.write("genomeDir=/shared/workspace/software/genomes/Hsapiens/star\n")
    filewriter.write("\n")
    filewriter.write("if [ ! -f $output_file/\"Aligned.out.sam\" ]; then\n")
    filewriter.write("   $STAR \\\n")
    filewriter.write("   --genomeDir $genomeDir \\\n")
    filewriter.write("   --readFilesIn $input_fastq \\\n")
    filewriter.write("   --outFileNamePrefix $output_file/ \\\n")

    for option in option_list:
        filewriter.write("   --" + option[0] + " " + option[1] + " \\\n")
    filewriter.write("fi\n")

if __name__ == "__main__":
    option_list = [["option_1", "my_option_1"], ["option_2", "my_option_2"]]
    configure_star(option_list)
