__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import paramiko
from scp import SCPClient

## connecting the master instance of a CFNCluster by the hostname, username and private key file
## return a ssh client
def connect_master(hostname, username, private_key_file):
    """Connects to a remote host.

    Uses a user's username, hostname, and private key file to 
    connect to username@hostname. Used to connect to cluster head node.

    Args:
        hostname: A string representing the IP address of the target host.
        username: A string representing the target username to connect to.
        private_key_file: A string path to the private key file needed
            to connect to the host.

    Returns:
        A paramiko SSHClient instance that stores the connection to
        username@hostname. <http://docs.paramiko.org/en/2.4/api/client.html>
    """
    private_key = paramiko.RSAKey.from_private_key_file(private_key_file)
    ssh_client = paramiko.SSHClient()
    ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    print("connecting")
    ssh_client.connect(hostname=hostname, username=username, pkey=private_key)
    print("connected")

    return ssh_client

## executing command using the ssh client
def execute_command(ssh_client, command,verbose=False):
    """Runs a command on a remote host.

    Used to execute scripts on the cluster from local controllers.
    Only the stderr and stdout from the scripts are recorded, 
    return values/exit statuses are not.

    Args:
        ssh_client: A paramiko SSHClient instance that describes
            the connection to the remote host.
        command: A string command to be run on the remote host.

    Returns:
        A string concatenation of the stdout and stderr
        of the command. Will not provide useful information
        if command does not output anything to stderr or
        stdout.
    """
    if verbose:
        print("Executing {}".format(command))
    stdin, stdout, stderr = ssh_client.exec_command(command)
    result = stdout.read().decode("utf-8")
    result += stderr.read().decode("utf-8")

    return result

def list_dir(ssh_client, directory):
    return execute_command(ssh_client, "ls {}".format(directory)).split()

def copy_file(ssh_client, localpath, remotepath):
    """Copies local file to remote host.

    Uses scp to copy a local file to specified location on 
    the remote host. Local file remains unchanged.

    Args:
        ssh_client: A paramiko SSHClient instance that describes
            the connection to the remote host.
        localpath: A string absolute path to the local file to be copied.
        remotepath: A string absolute path to the target location
            of the file to be copied.

    Returns:
        None
    """
    # SCPCLient takes a paramiko transport as its only argument
    scp = SCPClient(ssh_client.get_transport())
    print(localpath)
    print(remotepath)
    scp.put(localpath, remotepath)

def copy_gatk(ssh_client, localpath):
    """Copies GATK jar file to remote host.

    GATK licensing does allow for the distribution of its
    jar file in the cluster snapshot, so users must download
    it on their own and upload it to the cluster.

    Args:
        ssh_client: A paramiko SSHClient instance that describes
            the connection to the remote host.
        localpath: A string absolute path to the GATK jar file.

    Returns:
        None
    """
    # SCPCLient takes a paramiko transport as its only argument
    remotepath = "/shared/workspace/software/gatk/3.8-0/"
    copy_file(ssh_client, localpath, remotepath)

## close the ssh connection
def close_connection(ssh_client):
    """Closes connection to remote host.

    After this function the SSHClient instance
    is no longer able to perform anything on the remote host.
    If further operatinos are needed connect_master must
    be run again.

    Args:
        ssh_client: A paramiko SSHClient instance that describes
            the connection to the remote host.

    Returns:
        None
    """
    ssh_client.close()
