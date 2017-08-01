__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import yaml

## parseing yaml file
def parse_yaml_file(yaml_file):

    documents = {}
    samples = []
    stream = open(yaml_file, "r")
    docs = yaml.load_all(stream)
    for doc in docs:
        for k,v in doc.items():
            if k == "sample":
                samples.append(v)
            else:
                documents.update({k:v})
            print k, "->", v
        print "\n",
    documents.update({"sample":samples})

    return documents

if __name__ == "__main__":
    yaml_file = "/Users/guorongxu/Desktop/Sample_cDNA.yaml"
    documents = parse_yaml_file(yaml_file)

    print documents.get("samples")

