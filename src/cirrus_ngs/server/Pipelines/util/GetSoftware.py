import os

software_dir = "/shared/workspace/software/"
not_software = ["R-packages", "references", "new_homer", "kallisto_util", "ghostscript-9.19-linux-x86_64", "anaconda3"]

possible_tools = os.listdir(software_dir)

tools = [tool for tool in possible_tools if len(
    list(filter(lambda x : not x.startswith("."),os.listdir(software_dir + tool)))) == 1]

tools = [tool for tool in possible_tools if tool not in not_software]
tools_and_vers = {tool:list(filter(lambda x: not x.startswith("."), os.listdir(software_dir + tool))) for tool in tools}

result = {key:[v[v.index("-")+1:] + " (installed with {})".format(v[:v.index("-")]) if "pip" in v or "conda" in v
    else v for v in val] for key,val in tools_and_vers.items()}

#tools_and_vers = {key:[v.replace("anaconda-", "installed with conda ") for v in val] for key,val in tools_and_vers.items()}

print(result)
