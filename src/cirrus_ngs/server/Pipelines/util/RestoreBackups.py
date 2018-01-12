import os

config_dir = "/shared/workspace/Pipelines/config"
scripts_dir = "/shared/workspace/Pipelines/scripts"

for path,dirs,files in os.walk(config_dir):
    for curr_file in files:
        if curr_file.endswith("BACKUP"):
            os.rename("{}/{}".format(path, curr_file), "{}/{}".format(path, os.path.splitext(curr_file)[0]))

for path,dirs,files in os.walk(scripts_dir):
    dirs[:] = [d for d in dirs if not d == "deprecated"]
    for curr_file in files:
        if curr_file.endswith("BACKUP"):
            os.rename("{}/{}".format(path, curr_file), "{}/{}".format(path, os.path.splitext(curr_file)[0]))
