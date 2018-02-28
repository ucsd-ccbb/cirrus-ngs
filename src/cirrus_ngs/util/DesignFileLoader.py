__doc__ = """
Used in the pipeline jupyter notebooks to extract information
from the design file given by the user
"""
__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import re

## load design file
def load_design_file(design_file):
    """Parses the design file the user input into the jupyter notebook for the pipeline.
    
    See README.md for details on design file formats. This function should
    be called from a pipeline's jupyter notebook. 

    Args:
        design_file: string absolute path to the project's design file

    Returns:
        a tuple of form:
        (
            list of dictionaries containing information about each sample,
            list of all groups in the project,
            dictionary with key=normal/chip samples and values=tumor/input samples
        )
    """
    sample_list = []
    group_list = []
    normal_samples = {}
    tumor_samples = {}
    
    with open(design_file, 'r+') as f:
        for line in f:
            if not line.startswith("##"):
                fields = line.split("\t")

                #need at least two tab separated fields in design file
                if not len(fields) >= 2:
                    raise IndexError("""Design file lines must follow format
<sample_name><OPTIONAL sample_name_reverse_reads><TAB><group_name>.
No tabs were found in line:\n\t\"{}\" """.format(line.strip()))

                if fields[0].find(",") > -1:
                    sample_pair = fields[0].split(",")
                    paired_samples_1 = sample_pair[0]
                    paired_samples_2 = sample_pair[1]
                    sample_list.append([paired_samples_1, paired_samples_2])
                else:
                    sample_list.append([fields[0]])

                group_list.append(fields[1].rstrip())

                #if paired data is available
                if len(fields) == 3:
                    fields[2] = fields[2].rstrip()
                    if fields[2] == "Normal" or fields[2] == "Chip":
                        normal_samples[group_list[-1]] = sample_list[-1][0].split(".")[0]
                    elif fields[2] == "Tumor" or fields[2] == "Input":
                        tumor_samples[group_list[-1]] = sample_list[-1][0].split(".")[0]
                    else:
                        raise ValueError("Design file third column must be either \"Normal\"/\"Tumor\" or \"Chip\"/\"Input\". Current value is {}".format(fields[2]))

    pair_list = {normal_samples[x]: tumor_samples[x] for x in normal_samples}

    print(sample_list)
    print(group_list)
    print(pair_list)
    
    return sample_list, group_list, pair_list
