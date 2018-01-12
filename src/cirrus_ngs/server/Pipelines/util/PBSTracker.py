__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import re
import time
from subprocess import Popen, PIPE

def trackPBSQueue(minutes, shell_script):
    """Freezes calling function until all scripts with the same name are done running.

    This function waits until all of a given shell script
    has left the PBS queue before finishing. It is used to prevent
    moving forward in the pipeline before the prerequisite steps
    are completed. Uses qstat -j <shell_script> to check if
    the shell script is still running. After every check
    for completion it prints the number of jobs it's waiting for.

    Args:
        minutes: number of minutes to wait in between checks for completion
        shell_script: name of shell script to wait for, without .sh extension

    Returns:
        None
    """
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
