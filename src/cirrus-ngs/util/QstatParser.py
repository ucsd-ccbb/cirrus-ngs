__author__ = "Mustafa Guler<mguler@ucsd.edu>"

## parses qstat output to give user status information
def parse_qstat(qstat, job_name):
    #removes header line
    qstat[:] = qstat[2:]

    jobs = {"qw": 0,
            "r" : 0,
            "d" : 0,
            "E" : 0}

    for line in qstat:
        line = line.strip().split()
        if job_name in line[2]:
            if "E" in line[4]:
                jobs["E"] += 1
            elif "d" in line[4]:
                jobs["d"] += 1
            elif "t" in line[4]: #do nothing, catching edge
            else:
                jobs[line[4]] += 1

    print("The status of your \"%s\" job:" % job_name)
    print("\tThere are %d jobs currently running." % jobs["r"])
    print("\tThere are %d jobs currently queued." % jobs["qw"])
    if jobs["d"] > 0:
        print("\tThere are %d jobs that are hanging." % jobs["d"])
        print("\t\tPlease delete the hanging jobs and try again")
    if jobs["E"] > 0:
        print("\tThere are %d jobs that have encountered an error." % jobs["E"])
        print("\t\tTo see an error breakdown run the error checking cell.")
    
    zeros = [0] * len(jobs.values())
    if jobs.values() == zeros:
        print("Your \"%s\" job has finished!" % job_name)

