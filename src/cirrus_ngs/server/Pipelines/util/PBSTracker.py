__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import re
import time
from subprocess import Popen, PIPE

def trackPBSQueue(minutes, shell_script):

    qstat_pipe = Popen(["qstat", "-j", shell_script], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    qstat_out, qstat_err = qstat_pipe.communicate()

    while len(qstat_err.decode('utf-8')) == 0:

        is_done = True
        index = 0

        qstat_pipe = Popen(["qstat", "-j", shell_script], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        qstat_out, qstat_err = qstat_pipe.communicate()

        jobs = re.split('=+\n', qstat_out.decode('utf-8'))[1:]

        print("{} job(s) are running...".format(len(jobs)))
        time.sleep(minutes * 60)

    print()
    print("No jobs are running")
