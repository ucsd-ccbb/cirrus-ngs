import os
from collections import defaultdict

config_dir = "/shared/workspace/Pipelines/config/"

config_spec_dirs = filter(lambda x : os.path.isdir(config_dir + x), os.listdir(config_dir))

config_files = {pipeline:filter(lambda x : not x.startswith("."), os.listdir(config_dir + pipeline)) for pipeline in config_spec_dirs}
config_cats = defaultdict(dict)

for pipeline in config_files:
    for config_file in config_files[pipeline]:
        workflow = os.path.splitext("_".join(config_file.split("_")[1:]))[0]
        f = open(config_dir + "{}/{}".format(pipeline, config_file))
        config_cats[pipeline][workflow] = f.read()
        f.close()

print(dict(config_cats))
