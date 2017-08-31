__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import re
import time
from subprocess import Popen, PIPE

def trackPBSQueue(minutes, shell_script):
    #removes any directories from the start of shell_script and adds ext
    shell_script = shell_script.split("/")[-1] + ".sh"

    qstat_pipe = Popen(["qstat", "-j", shell_script], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    qstat_out, qstat_err = qstat_pipe.communicate()

    while True:
        time.sleep(minutes * 60)

        qstat_pipe = Popen(["qstat", "-j", shell_script], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        qstat_out, qstat_err = qstat_pipe.communicate()

        if not len(qstat_err.decode('utf-8')) == 0:
            break

        jobs = re.split('=+\n', qstat_out.decode('utf-8'))[1:]

        print("{} job(s) are running...".format(len(jobs)))

    print()
    print("No jobs are running")

if __name__ == "__main__":
    trackPBSQueue(1, "reallylongname.sh")
