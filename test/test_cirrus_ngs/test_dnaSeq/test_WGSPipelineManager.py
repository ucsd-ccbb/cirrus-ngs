import unittest
import sys
import os
sys.path.append(os.getcwd().replace("test", "src"))
sys.path.append(sys.path[-1] + "/cirrus_ngs")
import cirrus_ngs.dnaSeq.WGSPipelineManager as WGSPipelineManager
import cirrus_ngs.cfnCluster.ConnectionManager as ConnectionManager
import tempfile
import re

class test_WGSPipelineManager(unittest.TestCase):
    def test_check_status(self):
        hostname = "35.162.87.9"
        username = "ec2-user"
        key = "/home/mustafa/interns_oregon_key.pem"
        ssh_client = ConnectionManager.connect_master(hostname, username, key)

        temp = tempfile.NamedTemporaryFile(mode="w+", delete=False)
        temp.write("""#!/bin/bash
for i in `seq 1 10000`; do
    echo $(($i * $i))
done""")

        temp.close()

        localpath = temp.name
        remotepath = "/home/ec2-user"

        ConnectionManager.copy_file(ssh_client, localpath, remotepath)
        ConnectionManager.execute_command(ssh_client, "qsub -N tmp " + os.path.basename(localpath))
        job_details = WGSPipelineManager.check_status(ssh_client, "tmp")

        #checks to makes sure only one job was found
        self.assertEqual(len(job_details), 1)

        pattern = """=*\njob_number.*\nexec_file.*\nsubmission_time.*
owner.*\nuid.*\ngroup.*\ngid.*\nsge_o_home.*\nsge_o_log_name.*
sge_o_path.*\nsge_o_shell.*\nsge_o_workdir.*\nsge_o_host.*
account.*\nmail_list.*\nnotify.*\njob_name.*\njobshare.*
env_list.*\nscript_file.*\nbinding.*\njob_type.*\nscheduling info.*"""

        self.assertNotEqual(re.match(pattern, job_details[0]), None)
        ConnectionManager.execute_command(ssh_client, "rm tmp*")
        job_num = re.search("job_number:\s*(\d*)\n", job_details[0]).group(1)
        ConnectionManager.execute_command(ssh_client, "qdel %d" % int(job_num))
    
        
        


if __name__ == "__main__":
    unittest.main(module=__name__, buffer=True, exit=False)

