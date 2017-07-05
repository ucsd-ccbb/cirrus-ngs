__author__ = "Mustafa Guler<mguler@ucsd.edu>"

## parses qstat output to give user status information
def parse_qstat(qstat, job_name):
    #removes header line
    qstat[:] = qstat[2:]

    job_count = 0

    for line in qstat:
        line = line.strip().split()
        if job in line[2]:
            if line[4] == "qw":
                print("Your job " + job_name + " has been queued.")
            elif line[4] == "r":
                print("Your job " + job_name + " is currently running.")
                print("It was started at " + line[6] " on " + line[5])
            elif d in line[4]:
                print("Your job is hanging. Please delete these hanging jobs and run again.")
            elif E in line[4]:
                print("There has been an error with your job.")
            else
                print("If you can see this message I messed up.")








