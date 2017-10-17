from cfnCluster import ConnectionManager

def get_scripts(ssh_client):
    all_scripts = eval(ConnectionManager.execute_command(ssh_client, "python /shared/workspace/Pipelines/util/GetScripts.py"))
    #scripts = set(scripts.split("\n"))
    return all_scripts

def print_scripts(all_scripts):
    print("/shared/workspace/Pipelines/scripts:{}".format(all_scripts[0]))
    num_tabs = 1

    for pipeline in all_scripts[1]:
        print("\t"*num_tabs + "{}:".format(pipeline), end = "")

        num_tabs+=1
        curr_upper_scripts = [x for x in all_scripts[1][pipeline] if isinstance(x, str)]
        print(curr_upper_scripts)

        for workflow in all_scripts[1][pipeline]:
            if isinstance(workflow, dict):
                print("\t"*num_tabs,end="")
                print(workflow)
        num_tabs-=1

def print_target_scripts(all_scripts, pipeline, workflow):
    if not pipeline:
        print("scripts shared by all pipelines: {}".format(", ".join(all_scripts[0])))
        return
    if not workflow:
        print("scripts used by all {} workflows:".format(pipeline), end=" ")
        curr_upper_scripts = [x for x in all_scripts[1][pipeline] if isinstance(x, str)]
        print(", ".join(curr_upper_scripts))
        for workflow in all_scripts[1][pipeline]:
            if isinstance(workflow, dict):
                for key in workflow:
                    print("\t{}: {}".format(key, ", ".join(workflow[key])))
        return

    print("scripts under the {} {} workflow:".format(pipeline, workflow), end=" ")
    for entry in all_scripts[1][pipeline]:
        if isinstance(entry, dict) and entry.get(workflow, False):
            print(", ".join(list(map(str, entry[workflow]))))
