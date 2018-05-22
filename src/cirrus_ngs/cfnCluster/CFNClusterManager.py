__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import os
import re
from shutil import copyfile

## Global variable: default configuation file for CFNCluster
config_template_file = os.environ['HOME'] + "/.cfncluster/config" #os.getcwd().replace("notebooks", "data") + "/config"

def install_cfn_cluster():
    """Installs cfncluster framework.

    Uses pip to install the cfncluster framework and cli. 
    Does not perform updates or any changes if it
    is already installed with pip. Will not work if
    pip executable is not in $PATH. Prints output of 
    "pip install cfncluster" command to stdout.

    Args:
        None

    Returns:
        None
    """
    print("Installing cfncluster package...")
    print(os.popen("pip install cfncluster").read())

def upgrade_cfn_cluster():
    """Updates cfncluster framework.

    Uses pip to update cfncluster framework and cli. The
    cfncluster framework must already be installed on machine.
    Will not work if pip executable is not in $PATH. 
    Prints output of "pip install --upgrade cfncluster" command to 
    stdout.

    Args:
        None
    
    Returns:
        None
    """
    print("Upgrading cfncluster package...")
    print(os.popen("pip install --upgrade cfncluster").read())

def make_config_file():
    """Creates template cfncluster configuration file.

    Writes to the ~/.cfncluster/config file. If
    the configuraiton file directory does not exist
    it is created. If the configuration file does not 
    already exist it is created. 

    Args:
        None
    
    Returns:
        None
    """
    dir = os.path.dirname(config_template_file)
    if not os.path.exists(dir):
        os.makedirs(dir)

    filewriter = open(config_template_file, "w+")
    filewriter.write("[aws]" + "\n")

    filewriter.write("aws_region_name = us-east-1" + "\n")
    filewriter.write("aws_access_key_id = ***" + "\n")
    filewriter.write("aws_secret_access_key = ***" + "\n")

    filewriter.write("[cluster cfncluster]" + "\n")
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
    filewriter.write("#post_install = s3://bucket_name/path/to/postinstall.sh" + "\n")

    filewriter.write("[vpc ucsd]" + "\n")
    filewriter.write("master_subnet_id = subnet-00000000" + "\n")
    filewriter.write("vpc_id = vpc-00000000" + "\n")

    filewriter.write("[global]" + "\n")
    filewriter.write("update_check = true" + "\n")
    filewriter.write("sanity_check = true" + "\n")
    filewriter.write("cluster_template = cfncluster" + "\n")

    filewriter.write("[ebs custom]" + "\n")
    filewriter.write("ebs_snapshot_id = snap-a6e477ff" + "\n")
    filewriter.write("volume_size = 200" + "\n")
    filewriter.close()


## viewing CFNCluster configuration settings
def view_cfncluster_config():
    """Prints contents of cfncluster configuration file.

    If the file does not exist, it creates a template version
    of it in ~/.cfncluster/config using make_config_file. The contents are 
    printed to stdout.

    Args:
        None

    Returns:
        None
    """
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
            print(line[:-1])

## inserting AWS access keys
def insert_access_keys(aws_access_key_id="***",
                       aws_secret_access_key="***"):
    """Adds AWS key ID and keys to cfncluster configuration file.

    If the configuration file doesn't exist it creates the template 
    version first using make_config_file. 
    The "aws_access_key_id" field and "aws_secret_access_key"
    fields of the configuration file are rewritten based on 
    this functions args. 

    Args:
        aws_access_key_id: A string AWS access key ID.
        aws_secret_access_key: A string AWS secret access key.

    Returns:
        None
    """

    if not os.path.isfile(config_template_file):
        print("config_template_file not exist!")
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
    """Adds region to cfncluster configuration file.

    If the configuration file doesn't exist the template version is 
    created with make_config_file. The "aws_region_name" field of
    the configuration file is rewritten based on the function args. 

    Args:
        aws_region_name: A string valid AWS region. 

    Returns:
        None
    """
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
    """Adds private key name to cfncluster configuration file.

    If the configuration file doesn't exist the template version is
    created with make_config_file. The "key_name" field
    of the configuration file is rewritten based on the
    function args. 

    Args:
        key_name: A string absolute path to a private key file.

    Returns:
        None
    """
    if not os.path.isfile(config_template_file):
        make_config_file()

    #basename of private_key name, no .pem extension
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
    """Adds head node and computing node types to cfncluster configuration file.

    If the configuration file doesn't exist the template version is created 
    with make_config_file. The "master_instance_type" and "compute_instance_type"
    fields of the configuration file are rewritten based on the function args.

    Args:
        master_instance_type: A string valid AWS node type (e.g. "m3.large") used for
            the cluster's head node.
        compute_instance_type: A string valid AWS node type (e.g. "m3.large") used for
            the cluster's computing nodes. 

    Returns:
        None
    """
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
    """Adds starting cluster size to cfncluster configuration file.

    If the configuration file doesn't exist the template version is created 
    with make_config_file. The "initial_queue_size" field of the
    configuration file is rewritten based on the function args.

    Args:
        initial_cluster_size: A string initial number of compute nodes to 
            launch in cluster.

    Returns:
        None
    """
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
    """Adds spot price to cfncluster configuration file.

    If the configuration file doesn't exist the template version is created 
    with make_config_file. The "spot_price" field of the configuration file
    is rewritten based on the function args. 

    Args:
        spot_price: A string maximum spot price for the cluster.

    Returns:
        None
    """
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
    """Adds input and output s3 resources to cfncluster configuration file.

    If the configuration file doesn't exist the template version is created 
    with make_config_file. The "s3_read_resource" and "s3_read_write_resource"
    fields of the configuration file are rewritten based on the function
    args.

    Args:
        s3_read_resource: A string s3 bucket path that computing nodes 
            will have read access to.
        s3_read_write_resource: A string s3 bucket path that computing
            nodes will have read and write access to.

    Returns:
        None
    """
    if not os.path.isfile(config_template_file):
        make_config_file()
    ## s3://ucsd-ccbb-wgs-test-us-east-1/RNASeq_Pipeline_Code/test_data
    read_bucket_name = s3_read_resource[5:]
    read_bucket_name = read_bucket_name
    write_bucket_name = s3_read_write_resource[5:]
    write_bucket_name = write_bucket_name

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
    """Adds path to post installation script to cfncluster configuration file.

    If the configuration file doesn't exist the template version is created 
    with make_config_file. The "post_install" field of the configuration
    file is rewritten based on the function args. 

    Args:
        post_install: A string s3 bucket path to a script to be run
            on cluster after creation.

    Returns:
        None
    """
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
    """Add master subnet ID and VPC ID to the cfncluster configuration file.

    If the configuration file doesn't exist the template version is created 
    with make_config_file. The "master_subnet_id" and "vpc_id" fields of 
    the configuration file are rewritten based on the function args.

    Args:
        master_subnet_id: A string with format "subnet-{}" where {}
            is the ID of the subnet that the master server should
            be provisioned to.
        vpc_id: A string with format "vpc={}" where {}
            is the ID of the VPC the cluster should be
            provisioned to. 

    Returns:
        None
    """
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
    """Adds cirrus EBS snapshot ID to cfncluster configuration file

    If the configuration file doesn't exist the template version is created 
    with make_config_file. The "ebs_shapshot_id" field of the configuration
    file is rewritten based on the function args. The snapshot ID
    should correspond to a valid Cirrus snapshot.

    Args:
        ebs_snapshot_id: A string snapshot ID of a Cirrus cluster

    Returns:
        None
    """
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
    """Adds EBS volume size to the cfncluster configuration file.
    
    If the configuration file doesn't exist the template version is created 
    with make_config_file. The "volume_size" field of the configuration
    file is rewritten based on the function args. 

    Args:
        volume_size: A string representing the size of snapshot 
            volume in gigabytes

    Returns:
        None
    """
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
    """ Prints a list of current cfnclusters.

    This is based on the the cfncluster cli, it must be installed
    for the function to work. Prints the cluster names based
    on their configuration
    
    Args:
        None

    Returns:
        None
    """
    print(os.popen("cfncluster list").read())

## creating a CFNCluster with the specific cluster name
def create_cfn_cluster(cluster_name="mycluster"):
    """Creates a cfncluster.

    The cluster is created using the cfncluster
    cli, a configuration and a name. If a cluster
    with the given name already exists nothing is
    done.

    Args:
        cluster_name: A string name of the cluster to be created.

    Returns:
        A string IP address of the head node of the cluster. 
    """
    master_ip_address = ""

    response = os.popen("cfncluster status " + cluster_name).read()
    if response.find("CREATE_COMPLETE") > -1:
        print("cluster " + cluster_name + " does exist.")
        lines = re.split(r'\n+', response)
        for line in lines:
            if line.find("MasterPublicIP") > -1:
                master_ip_address = line[line.find("=") + 2:-1]
            print(line)

        return master_ip_address

    response = os.popen("cfncluster create " + cluster_name).read()

    lines = re.split(r'\n+', response)
    for line in lines:
        if line.find("MasterPublicIP") > -1:
            master_ip_address = line[line.find("=") + 2:-1]
        print(line)

    return master_ip_address

## deleting the specific CFNCluster
def delete_cfn_cluster(cluster_name="mycluster"):
    """Deletes a cfn cluster.

    This function completely gets rid of a cluster.
    Prints the output of the "cfncluster delete"
    call.

    Args:
        cluster_name: A string name of the cluster to be deleted
    
    Returns:
        None
    """
    print(os.popen("cfncluster delete " + cluster_name).read())

if __name__ == "__main__":
    view_cfncluster_config()
