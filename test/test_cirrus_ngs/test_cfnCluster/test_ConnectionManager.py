import unittest
import sys
import os
sys.path.append(os.getcwd().replace("test", "src"))
import cirrus_ngs.cfnCluster.ConnectionManager as ConnectionManager
import paramiko
import tempfile
import re

class test_ConnectionManager(unittest.TestCase):
    def test_paramiko(self):
        key_file = tempfile.NamedTemporaryFile()
        key_file.write(b"notakey")

        self.assertRaises(paramiko.SSHException, paramiko.RSAKey.from_private_key_file, key_file.name)

        key_file.close()
        new_key = "/home/mustafa/interns_oregon_key.pem"

        #checks to make sure a real key file works. will not be portable
        #leaving my ssh key for users to download for tests seems not smart
        paramiko.RSAKey.from_private_key_file(new_key)

    def test_connect_master(self):
        hostname = "35.162.87.9"
        username = "ec2-user"
        key_file = tempfile.NamedTemporaryFile()
        key_file.write(b"not_a_key")
        key_file.seek(0)

        self.assertRaises(paramiko.SSHException, ConnectionManager.connect_master, hostname, username, key_file.name)
        key_file.close()

        #this won't even work elsewhere but I don't want to put my keyfile into the eepo
        new_key = "/home/mustafa/interns_oregon_key.pem"
        ConnectionManager.connect_master(hostname, username, new_key)

        #checks if last line in the standard output is "connected"
        out = sys.stdout.getvalue().strip()
        last_line = out.split()[-1]
        self.assertEqual(last_line, "connected")

        #checks that connected and connecting only are printed once exactly
        num_connected = len(re.findall("connected", out))
        self.assertEqual(1, num_connected)

        num_connecting = len(re.findall("connecting", out))
        self.assertEqual(1, num_connecting)

    def test_execute_command(self):
        hostname = "35.162.87.9"
        username = "ec2-user"
        key = "/home/mustafa/interns_oregon_key.pem"
        ssh_client = ConnectionManager.connect_master(hostname, username, key)
        command = "pwd"

        #checks that the pwd command worked
        self.assertEqual(ConnectionManager.execute_command(ssh_client, command), "/home/ec2-user\n")

        ssh_client = "not an ssh_client"

        #makes sure that an error is raised when a non sshclient is passed in
        self.assertRaises(AttributeError, ConnectionManager.execute_command, ssh_client, command)

    def test_copy_file(self):
        hostname = "35.162.87.9"
        username = "ec2-user"
        key = "/home/mustafa/interns_oregon_key.pem"
        ssh_client = ConnectionManager.connect_master(hostname, username, key)
        
        temp = tempfile.NamedTemporaryFile()
        localpath = temp.name 
        remotepath = "/home/ec2-user"
        ConnectionManager.copy_file(ssh_client, localpath, remotepath)
        out = sys.stdout.getvalue().strip().split()[-2:]

        #checks that the copy file prints the local and remote paths
        self.assertEqual(out, [localpath, remotepath])

        ls_output = ConnectionManager.execute_command(ssh_client, 
                "ls tmp* | wc -l")
        ConnectionManager.execute_command(ssh_client, "rm tmp*")

        #checks that there is exactly 1 tempfile in the home directory of the server
        self.assertEqual(ls_output.strip(), "1")

        #makes sure it doesn't work with a nonfile
        self.assertRaises(FileNotFoundError, ConnectionManager.copy_file, 
                ssh_client, "fakefile", "/home/ec2-user")
    
    #########################################################################
    #copy_gatk, list_dir, and close_connection are considered trivial methods 
    #and are not tested
    #########################################################################

if __name__ == "__main__":
    unittest.main(module=__name__, buffer=True, exit=False)
