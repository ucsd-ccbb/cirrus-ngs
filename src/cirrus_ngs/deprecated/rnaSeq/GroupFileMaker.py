__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

## make a group file for analysis of Small RNA Sequencing pipeline
def make_group_file(workflow, project_name, group_file, sample_list):
    localpath = "/shared/workspace/data_archive/RNASeq/" + project_name + "/" + workflow + "/"

    group_table = {}
    for sample in sample_list:
        group = sample.get("group")
        if group not in group_table:
            group_table.update({group: 1})
        else:
            number = group_table.get(group)
            group_table.update({group: number + 1})

    filewriter = open(group_file, "w")

    filewriter.write(localpath + "all_gene_counts.txt\t3\n")
    for group in group_table:
        filewriter.write(group + "\t" + str(group_table.get(group)) + "\n")

    filewriter.close()
