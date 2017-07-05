__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import paramiko
from scp import SCPClient

## connecting the master instance of a CFNCluster by the hostname, username and private key file
## return a ssh client
def connect_master(hostname, username, private_key_file):
    private_key = paramiko.RSAKey.from_private_key_file(private_key_file)
    ssh_client = paramiko.SSHClient()
    ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    print "connecting"
    ssh_client.connect(hostname=hostname, username=username, pkey=private_key)
    print "connected"

    return ssh_client

## executing command using the ssh client
def execute_command(ssh_client, command):
    print "Executing {}".format(command)
    stdin , stdout, stderr = ssh_client.exec_command(command)
    print stdout.read()
    print stderr.read()

def copy_file(ssh_client, localpath, remotepath):
    # SCPCLient takes a paramiko transport as its only argument
    scp = SCPClient(ssh_client.get_transport())
    print localpath
    print remotepath
    scp.put(localpath, remotepath)

def copy_gatk(ssh_client, localpath):
    # SCPCLient takes a paramiko transport as its only argument
    scp = SCPClient(ssh_client.get_transport())
    remotepath = "/shared/workspace/software/gatk/GenomeAnalysisTK-3.3-0"
    print "coping " + localpath + " to cluster..."
    scp.put(localpath, remotepath)

## close the ssh connection
def close_connection(ssh_client):
    ssh_client.close()
