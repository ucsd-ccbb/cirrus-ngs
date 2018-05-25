import os
import sys
import shutil
import pkg_resources
import yaml
import string
import json
import re
import subprocess

def _clean_root_path(root_path):
    root_path = os.path.abspath(root_path)
    root_path = os.path.expanduser(root_path)
    return root_path

def _check_curr_dir():
    cirrus_config = "{}/.cirrus.yaml".format(os.getcwd())
    if os.path.exists(cirrus_config):
        return cirrus_config
    else:
        raise FileNotFoundError("{} is not a cirrus-ngs directory".format(os.getcwd()))

def init(root_path):
    root_path = _clean_root_path(root_path) 

    if not os.path.exists(root_path):
        os.makedirs(root_path)

    notebook_dir = "{}/notebooks".format(root_path)
    design_dir = "{}/design_files".format(root_path)
    cirrus_config = "{}/.cirrus.yaml".format(root_path)
    
    if not os.path.exists(notebook_dir):
        os.mkdir(notebook_dir)
    if not os.path.exists(design_dir):
        os.mkdir(design_dir)
    if not os.path.exists(cirrus_config):
        config_template = pkg_resources.resource_filename("cirrusngs", "cirrus.yaml")
        shutil.copy(config_template, cirrus_config)

    # copy cfn cluster setup notebook
    if not os.path.exists(os.path.join(notebook_dir, 
        "BasicCFNClusterSetup.ipynb")):
        notebook_templates = pkg_resources.resource_filename("cirrusngs", "notebooks")
        shutil.copy(os.path.join(notebook_templates, "BasicCFNClusterSetup.ipynb"),
                os.path.join(notebook_dir, "BasicCFNClusterSetup.ipynb"))

def add(pipeline_name):
    _check_curr_dir()

    # preprocess pipeline name
    pipeline_name = pipeline_name.lower()
    table = str.maketrans(dict.fromkeys(string.punctuation))
    pipeline_name = pipeline_name.translate(table)
    if not (pipeline_name.endswith("seq") or pipeline_name in {"wgs", "wes", "all"}):
        pipeline_name += "seq"

    # get possible pipeline names
    notebook_dir = pkg_resources.resource_filename("cirrusngs", "notebooks")
    design_dir = pkg_resources.resource_filename("cirrusngs", "design_files")
    notebooks = os.listdir(notebook_dir)
    notebooks = list(filter(lambda x:x.endswith(".ipynb"), notebooks))
    possible_names = list(map(lambda x:x.split("Template")[0].lower(), notebooks))

    # plain copying of all notebooks and design files, can't config them later
    if pipeline_name == "all":
        for nb in notebooks:
            if not nb.endswith("SeqTemplate.ipynb"):
                continue
            if nb.startswith("DNA"):
                target = "WGS|WESPipeline.ipynb"
            else:
                target = nb.replace("Template", "Pipeline")

            shutil.copy(os.path.join(notebook_dir, nb),
                    "notebooks/{}".format(target))
            target = nb.replace("Template.ipynb", "_design_example.txt")
            shutil.copy(os.path.join(design_dir, target),
                    "design_files/{}".format(target))
        return


    if not pipeline_name in possible_names:
        raise ValueError("{} is not a supported pipeline".format(pipeline_name))

    # copy files over
    target = notebooks[possible_names.index(pipeline_name)]
    shutil.copy("{}/{}".format(notebook_dir, target), 
            "notebooks/{}".format(target))
    target = target.replace("Template.ipynb", "_design_example.txt")
    shutil.copy("{}/{}".format(design_dir, target), 
            "design_files/{}".format(target))

def make_config(config_dict, config_file):
    _check_curr_dir()
    with open(config_file, "w") as f:
        f.write(yaml.dump(config_dict))

def apply_config(config_file, notebooks):
    _check_curr_dir()

    # possible values from user
    possible_nbs = os.listdir("notebooks")
    possible_nbs = list(filter(lambda x:x.endswith("Template.ipynb"), possible_nbs))
    possible_names = list(map(lambda x:x.split("Template")[0].lower(), possible_nbs))

    # load config file they want to apply
    with open(config_file) as f:
        config_dict = yaml.load(f)

    for notebook in notebooks:
        # normalize their notebook name
        notebook = notebook.lower()
        table = str.maketrans(dict.fromkeys(string.punctuation))
        notebook = notebook.translate(table)
        if not (notebook.endswith("seq") or notebook == "WGS" or notebook == "WES"):
            notebook += "seq"

        if not notebook in possible_names:
            raise ValueError("{} is not a template in your cirrus notebooks directory"
                    .format(notebook))

        target = possible_nbs[possible_names.index(notebook)]

        new_nb = os.path.join("notebooks", "{}_{}".format(
            config_dict["project_name"], target.replace("Template","Pipeline")))
        old_nb = os.path.join("notebooks", target)

        # load old notebook
        with open(old_nb) as f:
            curr_nb = json.load(f)

        # find first cell with code in it
        code_ind = 0
        for cell in curr_nb["cells"]:
            if cell["cell_type"] == "code":
                break
            code_ind += 1

        # format first code cell with config file
        source = curr_nb["cells"][code_ind]["source"]
        source = list(map(lambda x:
            x.format_map(config_dict) if '{' in x and '}' in x else x, source))
        curr_nb["cells"][code_ind]["source"] = source

        # write new notebook
        with open(new_nb, "w") as f:
            json.dump(curr_nb, f)
