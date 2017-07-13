import unittest
import sys
import os
sys.path.append(os.getcwd().replace("test", "src"))
sys.path.append(sys.path[-1] + "/cirrus_ngs")
import cirrus_ngs.dnaSeq.WGSPipelineManager as WGSPipelineManager
import cirrus_ngs.cfnCluster.ConnectionManager as ConnectionManager
import cirrus_ngs.util.QstatParser as QstatParser
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
        WGSPipelineManager.check_status(ssh_client, "tmp")
        job_result = sys.stdout.getvalue().strip().split("\n")
        for i in range(len(job_result)):
            if job_result[i] == "":
                job_result = job_result[i+1:]
                break

        correct_result = """The status of your "tmp" job:
There are 0 jobs currently running.
There are 1 jobs currently queued.
\tQueued jobs:
\t\ttmp"""

        qstat = ConnectionManager.execute_command(ssh_client, "qstat")
        job_num = QstatParser.get_job_ids(qstat)[0][0]
        ConnectionManager.execute_command(ssh_client, "qdel {}".format(int(job_num)))
        ConnectionManager.execute_command(ssh_client, "rm tmp*")
    

        #checks to makes sure stdout matches intention
        self.assertEqual("\n".join(job_result), correct_result)

    def test_execute(self):
        self.fail()


if __name__ == "__main__":
    unittest.main(module=__name__, buffer=True, exit=False)

