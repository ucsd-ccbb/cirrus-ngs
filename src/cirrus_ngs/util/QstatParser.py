__author__ = "Mustafa Guler<mguler@ucsd.edu>"
import re
import os
import sys

def get_job_ids(qstat):
    if qstat == "\n":
        return None

    qstat_list = qstat.split("\n")

    #removes header lines and last entry from the extra newline
    qstat_list = qstat_list[2:-1]

    job_ids = [[int(qstat_list[x].split()[0]), qstat_list[x].split()[4]] 
            for x in range(len(qstat_list))]
    return job_ids

## parses qstat output to give user status information
def parse_qstat(job_info , job_name, logs):

    jobs = {"qw": [],
            "r" : [],
            "d" : [],
            "E" : []}

    for job in job_info:
        curr_job_name = re.search("job_name:\s*(\S*)", job[2]).group(1)
        if job_name in curr_job_name:
            if "E" in job[1]:
                jobs["E"].append(curr_job_name)
            elif "d" in job[1]:
                jobs["d"].append(curr_job_name)
            elif job[1] == "qw" or job[1] == "r":
                jobs[job[1]].append(curr_job_name)

    print("\nThe status of your \"{}\" job:".format(job_name))

    print("There are {} jobs currently running.".format(len(jobs["r"])))
    if len(jobs["r"]) > 0:
        print("\tRunning jobs:")
    for job in jobs["r"]:
        print("\t\t{}".format(job))
    print("There are {} jobs currently queued.".format(len(jobs["qw"])))
    if len(jobs["qw"]) > 0:
        print("\tQueued jobs:")
    for job in jobs["qw"]:
        print("\t\t{}".format(job))

    if len(jobs["d"]) > 0:
        print("\tThere are {} jobs that are hanging.".format(len(jobs["d"])))
        print("\t\tPlease delete the hanging jobs and try again")
    if len(jobs["E"]) > 0:
        print("\tThere are {} jobs that have encountered an error.".format(len(jobs["E"])))
        #print("\t\tTo see an error breakdown run the error checking cell.")

    job_done = False
    for log in logs:
        if job_name in log:
            job_done = True
            break

    zeros = {"qw": [],
             "r" : [],
             "d" : [],
             "E" : []}
    if jobs == zeros and job_done:
        print("Your \"{}\" job has finished!".format(job_name))
    elif jobs == zeros:
        print("Your \"{}\" job has not started yet.".format(job_name))

