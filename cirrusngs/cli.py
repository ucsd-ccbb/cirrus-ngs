import cirrusngs
import argparse
import yaml
import os
import sys

cirrus_desc = ("cirrus-ngs is a cloud-optimized primary analysis pipeline for "
"RNA-seq, miRNA-seq, ChIP-seq, and variant calling in whole-genome/whole-exome DNA-seq")

def init(args):
    cirrusngs.init(args.root_path)
    sys.stderr.write("cirrus-ngs has been initialized in {}\n".format(args.root_path))
    sys.stderr.write("""add notebooks and design file templates for different pipelines with:
    cirrus-ngs add <pipeline>\n""")
    

def config(args):
    if args.template and args.apply:
        print(args.template)
        print()
        raise argparse.ArgumentTypeError("Cannot use the '-t' and '-a' flags together in config command")

    cirrus_config = "{}/.cirrus.yaml".format(os.getcwd())
    template = cirrus_config
    if args.template:
        template = args.template

    with open(template) as f:
        config_dict = yaml.load(f)

    if not args.apply:
        # asks user for config values
        for var in config_dict:
            print("{} [{}]: ".format(var, config_dict[var]), end="")
            new_val = input()
            if new_val:
                config_dict[var] = new_val
        cirrusngs.make_config(config_dict, args.file_name)
        sys.stderr.write("Configuration has been saved to {}\n".format(args.file_name))
    else:
        cirrusngs.apply_config(args.file_name, args.apply)

def add(args):
    for pipeline_name in args.pipeline_names:
        cirrusngs.add(pipeline_name)

def main():
    parser = argparse.ArgumentParser(description=cirrus_desc, prog="cirrus-ngs")
    subs = parser.add_subparsers()
    subs.required = True
    subs.dest = "command"
    
    init_parser = subs.add_parser("init", help="initialize a cirrus project",
            description="initializes target directory to be a cirrus project directory")
    init_parser.add_argument("root_path", type=str, nargs="?", 
            default=os.path.abspath("."), help="root directory of cirrus project, defaults to current directory")
    init_parser.set_defaults(func=init)

    config_parser = subs.add_parser("config", 
            help="configure notebook fields, can create new config files or apply existing configs to notebook templates",
            description=("create configuration files, can supply templates with the '-t' option. "
                "apply configuration files to notebook templates using '-a' and a list of notebook "
                "template files. DO NOT use '-t' and '-a' at the same time"))
    config_parser.add_argument("-t", "--template", 
            help="template to base config file on, do not use with '-a' option")
    config_parser.add_argument("-a", "--apply", nargs="+",
            help="notebooks to apply this config file to, do not use with '-t' option")
    config_parser.add_argument("file_name", type=str,
            help="name of config file")
    config_parser.set_defaults(func=config)

    add_parser = subs.add_parser("add", help="add notebook and design file for a pipeline",
            description=("adds notebook and design file templates to cirrus project in current directory"
                "the pipeline_names should not be file names: use 'rna', 'chip', 'mirna', or 'dna' instead."
                "capitilization and adding 'seq' to the end of these does not make a difference."
                "'wgs' and 'wes' can be used instead of 'dna'."
                "alternatively, 'all' can be used if simplicity is desired. this will copy all notebooks over"))
    add_parser.add_argument("pipeline_names", nargs="+", type=str, 
            help="pipeline names to add files for, punctuation and caps don't matter (Example: chip rna WGS)")
    add_parser.set_defaults(func=add)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
