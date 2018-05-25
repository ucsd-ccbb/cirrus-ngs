import os
import sys

scripts_path = "/shared/workspace/Pipelines/scripts"
banned_dirs = ("deprecated", "docs", "old_scripts")


def get_scripts(scripts_path, banned_dirs):
    result = {}

    for root,dirs,files in os.walk(scripts_path):
        dirs[:] = [d for d in dirs if not os.path.basename(d) in banned_dirs]
        curr_dir = os.path.basename(root)
        inner_dirs = list(map(os.path.basename, dirs))
        curr_scripts = list(map(os.path.basename, files))

        if curr_dir == "scripts":
            result.update({"All Pipelines":curr_scripts})
            result.update({pipeline:{} for pipeline in inner_dirs})
            continue

        if curr_dir in result:
            result[curr_dir] = {workflow:[] for workflow in inner_dirs}
            result[curr_dir].update({"All Workflows":files})
        else:
            curr_pipeline = list(filter(lambda x : curr_dir in result[x], result))[0]
            result[curr_pipeline][curr_dir] = curr_scripts
        
    return result
       
if __name__ == "__main__":
    if len(sys.argv) >= 3:
        print(get_scripts(sys.argv[1], sys.argv[2]))
    else:
        print(get_scripts(scripts_path, banned_dirs))

