import yaml

yaml_file = open("tools.yaml", "r")

doc = yaml.load(yaml_file)
print(doc)
