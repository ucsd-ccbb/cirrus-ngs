__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import re
import time
import subprocess

def trackPBSQueue(minutes, shell_script):

    jobs = 0
    while True:
        time.sleep(minutes * 60)

        is_done = True
        index = 0

        cmd = 'qstat'
        stdout = subprocess.check_output([cmd])

        lines = re.split(r'\n+', stdout)

        for line in lines:
            if line.find(shell_script) > -1 and not isDone(line):
                is_done = False
                index = index + 1

        if index > 0:
            if index != jobs:
                jobs = index
                print ("")
                print (str(index) , " job(s) are running...")

        if is_done:
            print ("")
            print ("No jobs are running...")
            break

def isDone(qstat_job):
    status_list = qstat_job.strip().split()
    if len(status_list) == 1:
        return True

    qstat_job_status = status_list[4]

    if re.match("d.*", qstat_job_status):
        return True
    return False

if __name__ == "__main__":
    trackPBSQueue(1, "PBS_Tracker")
