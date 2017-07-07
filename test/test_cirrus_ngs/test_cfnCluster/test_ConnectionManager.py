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

        #this is going to have to not be included on release. for obvious reasons
        new_key = "/home/mustafa/interns_oregon_key.pem"
        ConnectionManager.connect_master(hostname, username, new_key)

        #checks if last line in the standard output is "connected"
        out = sys.stdout.getvalue().strip()
        last_line = out.split()[-1]
        self.assertEqual(last_line, "connected")

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

        self.assertEqual(ConnectionManager.execute_command(ssh_client, command), "/home/ec2-user\n")

        ssh_client = "not an ssh_client"

        self.assertRaises(AttributeError, ConnectionManager.execute_command, ssh_client, command)

if __name__ == "__main__":
    unittest.main(module=__name__, buffer=True, exit=False)
