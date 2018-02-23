import os

software_dir = "/shared/workspace/software/"
not_software = ["R-packages", "references", "kallisto_util", "anaconda3"]

possible_tools = os.listdir(software_dir)

tools = [tool for tool in possible_tools if len(
    list(filter(lambda x : not x.startswith("."),os.listdir(software_dir + tool)))) == 1]

tools = [tool for tool in possible_tools if tool not in not_software]
tools_and_vers = {tool:list(filter(lambda x: not x.startswith("."), os.listdir(software_dir + tool))) for tool in tools}

result = {key:[v[v.index("-")+1:] + " (installed with {})".format(v[:v.index("-")]) if "pip" in v or "conda" in v
    else v for v in val] for key,val in tools_and_vers.items()}

print(result)
