__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import sys
import YamlFileReader


## make a group file for analysis of the pipeline
# group_file: path to group file to be generated
# yaml_file: path to yaml file for the project
# s3_path: directory for the all_gene_counts file on s3
def make_group_file(group_file, yaml_file, s3_path):

    documents = YamlFileReader.parse_yaml_file(yaml_file)
    sample_list = documents.get("sample")

    group_table = {}
    for sample in sample_list:
        group = sample.get("group")
        if group not in group_table:
            group_table.update({group: 1})
        else:
            number = group_table.get(group)
            group_table.update({group: number + 1})

    filewriter = open(group_file, "w")

    filewriter.write(s3_path + "/all_gene_counts.txt\t3\n")
    for group in group_table:
        filewriter.write(group + "\t" + str(group_table.get(group)) + "\n")

    filewriter.close()


if __name__ == "__main__":

    group_file = sys.argv[1]
    yaml_file = sys.argv[2]
    sample_list = sys.argv[3]

    make_group_file(group_file, yaml_file, sample_list)

