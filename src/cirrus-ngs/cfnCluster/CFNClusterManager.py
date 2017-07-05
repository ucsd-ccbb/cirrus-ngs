__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import os
import re
from shutil import copyfile

## Global variable: default configuation file for CFNCluster
config_template_file = os.getcwd().replace("notebooks", "data") + "/config"

def install_cfn_cluster():
    print "Installing cfncluster package..."
    print os.popen("pip install cfncluster").read()

def upgrade_cfn_cluster():
    print "Upgrading cfncluster package..."
    print os.popen("pip install --upgrade cfncluster").read()

def make_config_file():
    dir = os.path.dirname(config_template_file)
    if not os.path.exists(dir):
        os.makedirs(dir)

    filewriter = open(config_template_file, "w+")
    filewriter.write("[aws]" + "\n")

    filewriter.write("aws_region_name = us-east-1" + "\n")
    filewriter.write("aws_access_key_id = ***" + "\n")
    filewriter.write("aws_secret_access_key = ***" + "\n")

    filewriter.write("[cluster elasticsearch]" + "\n")
    filewriter.write("vpc_settings = ucsd" + "\n")
    filewriter.write("key_name = ***" + "\n")
    filewriter.write("master_instance_type = m3.large" + "\n")
    filewriter.write("compute_instance_type = r3.4xlarge" + "\n")
    filewriter.write("initial_queue_size = 0" + "\n")
    filewriter.write("cluster_type = spot" + "\n")
    filewriter.write("spot_price = 0.5" + "\n")
    filewriter.write("ebs_settings = custom" + "\n")

    filewriter.write("s3_read_resource = arn:aws:s3:::bucket_name" + "\n")
    filewriter.write("s3_read_write_resource = arn:aws:s3:::bucket_name/*" + "\n")
    filewriter.write("post_install = s3://bucket_name/path/to/postinstall.sh" + "\n")

    filewriter.write("[vpc ucsd]" + "\n")
    filewriter.write("master_subnet_id = subnet-00000000" + "\n")
    filewriter.write("vpc_id = vpc-00000000" + "\n")

    filewriter.write("[global]" + "\n")
    filewriter.write("update_check = true" + "\n")
    filewriter.write("sanity_check = true" + "\n")
    filewriter.write("cluster_template = elasticsearch" + "\n")

    filewriter.write("[ebs custom]" + "\n")
    filewriter.write("ebs_snapshot_id = snap-a6e477ff" + "\n")
    filewriter.write("volume_size = 200" + "\n")
    filewriter.close()

## viewing CFNCluster configuration settings
def view_cfncluster_config():
    if not os.path.isfile(config_template_file):
        make_config_file()

    cfncluster_config = os.environ['HOME'] + "/.cfncluster/config"

    if not os.path.isfile(cfncluster_config):
        if not os.path.exists(os.path.dirname(cfncluster_config)):
            os.makedirs(os.path.dirname(cfncluster_config))
        copyfile(config_template_file, cfncluster_config)

    with open(config_template_file) as fp:
        lines = fp.readlines()
        for line in lines:
            print line[:-1]

## inserting AWS access keys
def insert_access_keys(aws_access_key_id="***",
                       aws_secret_access_key="***"):

    if not os.path.isfile(config_template_file):
        make_config_file()

    with open(config_template_file, 'r+') as f:
        lines = f.readlines()
        f.seek(0)
        f.truncate()
        for line in lines:
            if 'aws_access_key_id' in line:
                line = line.replace(line[line.find("=") + 2:], aws_access_key_id + "\n")

            if 'aws_secret_access_key' in line:
                line = line.replace(line[line.find("=") + 2:], aws_secret_access_key + "\n")

            f.write(line)

## configuring aws region name
def config_aws_region_name(aws_region_name="us-east-1"):
    if not os.path.isfile(config_template_file):
        make_config_file()

    with open(config_template_file, 'r+') as f:
        lines = f.readlines()
        f.seek(0)
        f.truncate()
        for line in lines:
            if 'aws_region_name' in line:
                line = line.replace(line[line.find("=") + 2:], aws_region_name + "\n")

            f.write(line)

## configuring key pem file
def config_key_name(key_name):
    if not os.path.isfile(config_template_file):
        make_config_file()

    private_key = key_name[key_name.rfind("/") + 1:-4]

    with open(config_template_file, 'r+') as f:
        lines = f.readlines()
        f.seek(0)
        f.truncate()
        for line in lines:
            if 'key_name' in line:
                line = line.replace(line[line.find("=") + 2:], private_key + "\n")

            f.write(line)

## configuring master instance type and computer instance types
def config_instance_types(master_instance_type="m3.large", compute_instance_type="r3.2xlarge"):
    if not os.path.isfile(config_template_file):
        make_config_file()

    with open(config_template_file, 'r+') as f:
        lines = f.readlines()
        f.seek(0)
        f.truncate()
        for line in lines:
            if 'master_instance_type' in line:
                line = line.replace(line[line.find("=") + 2:], master_instance_type + "\n")

            if 'compute_instance_type' in line:
                line = line.replace(line[line.find("=") + 2:], compute_instance_type + "\n")

            f.write(line)

## configuring initial cluster size
def config_initial_cluster_size(initial_cluster_size="1"):
    if not os.path.isfile(config_template_file):
        make_config_file()

    with open(config_template_file, 'r+') as f:
        lines = f.readlines()
        f.seek(0)
        f.truncate()
        for line in lines:
            if 'initial_queue_size' in line:
                line = line.replace(line[line.find("=") + 2:], initial_cluster_size + "\n")

            f.write(line)

## configuring spot price for computer instances
def config_spot_price(spot_price="0.5"):
    if not os.path.isfile(config_template_file):
        make_config_file()

    with open(config_template_file, 'r+') as f:
        lines = f.readlines()
        f.seek(0)
        f.truncate()
        for line in lines:
            if 'spot_price' in line:
                line = line.replace(line[line.find("=") + 2:], spot_price + "\n")

            f.write(line)

## configuring S3 read/write resource: bucket name
def config_s3_resource(s3_read_resource="s3://bucket_name/", s3_read_write_resource="s3://bucket_name/"):
    if not os.path.isfile(config_template_file):
        make_config_file()
    ## s3://ucsd-ccbb-wgs-test-us-east-1/RNASeq_Pipeline_Code/test_data
    read_bucket_name = s3_read_resource[5:]
    read_bucket_name = read_bucket_name[:read_bucket_name.find("/")]
    write_bucket_name = s3_read_write_resource[5:]
    write_bucket_name = write_bucket_name[:write_bucket_name.find("/")]

    with open(config_template_file, 'r+') as f:
        lines = f.readlines()
        f.seek(0)
        f.truncate()
        for line in lines:
            if 's3_read_resource' in line:
                line = line.replace(line[line.find("=") + 2:], "arn:aws:s3:::" + read_bucket_name + "\n")

            if 's3_read_write_resource' in line:
                line = line.replace(line[line.find("=") + 2:], "arn:aws:s3:::" + write_bucket_name + "/*\n")

            f.write(line)

## configuring post installation shell script for creating CFNCluster
def config_post_install(post_install="s3://bucket_name/path/to/postinstall.sh"):
    if not os.path.isfile(config_template_file):
        make_config_file()

    with open(config_template_file, 'r+') as f:
        lines = f.readlines()
        f.seek(0)
        f.truncate()
        for line in lines:
            if 'post_install' in line:
                line = line.replace(line[line.find("=") + 2:], post_install + "\n")

            f.write(line)

## configuring vpc and subnet ids
def config_vpc_subnet_id(master_subnet_id="subnet-00000000",
                       vpc_id="vpc-00000000"):
    if not os.path.isfile(config_template_file):
        make_config_file()

    with open(config_template_file, 'r+') as f:
        lines = f.readlines()
        f.seek(0)
        f.truncate()
        for line in lines:
            if 'master_subnet_id' in line:
                line = line.replace(line[line.find("=") + 2:], master_subnet_id + "\n")

            if 'vpc_id' in line:
                line = line.replace(line[line.find("=") + 2:], vpc_id + "\n")

            f.write(line)

## configuring EBS snapshot id
def config_ebs_snapshot_id(ebs_snapshot_id="snap-a6e477ff"):
    if not os.path.isfile(config_template_file):
        make_config_file()

    with open(config_template_file, 'r+') as f:
        lines = f.readlines()
        f.seek(0)
        f.truncate()
        for line in lines:
            if 'ebs_snapshot_id' in line:
                line = line.replace(line[line.find("=") + 2:], ebs_snapshot_id + "\n")

            f.write(line)

## configuring EBS volume size to attach to CFNCluster
def config_volume_size(volume_size="200"):
    if not os.path.isfile(config_template_file):
        make_config_file()

    with open(config_template_file, 'r+') as f:
        lines = f.readlines()
        f.seek(0)
        f.truncate()
        for line in lines:
            if 'volume_size' in line:
                line = line.replace(line[line.find("=") + 2:], volume_size + "\n")

            f.write(line)

## listing all current cfnclusters
def list_cfn_cluster():
    print os.popen("cfncluster list").read()

## creating a CFNCluster with the specific cluster name
def create_cfn_cluster(cluster_name="mycluster"):
    master_ip_address = ""

    response = os.popen("cfncluster status " + cluster_name).read()
    if response.find("CREATE_COMPLETE") > -1:
        print "cluster " + cluster_name + " does exist."
        lines = re.split(r'\n+', response)
        for line in lines:
            if line.find("MasterPublicIP") > -1:
                master_ip_address = line[line.find("=") + 2:-1]
            print line

        return master_ip_address

    response = os.popen("cfncluster create " + cluster_name).read()

    lines = re.split(r'\n+', response)
    for line in lines:
        if line.find("MasterPublicIP") > -1:
            master_ip_address = line[line.find("=") + 2:-1]
        print line

    return master_ip_address

## deleting the specific CFNCluster
def delete_cfn_cluster(cluster_name="mycluster"):
    print os.popen("cfncluster delete " + cluster_name).read()

if __name__ == "__main__":
    view_cfncluster_config()
