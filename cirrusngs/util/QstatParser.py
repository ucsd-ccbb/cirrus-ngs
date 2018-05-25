import yaml
from collections import defaultdict
from cirrusngs.managers import ConnectionManager
import datetime
import re

def check_status(ssh_client, step_name, pipeline, workflow, project_name,analysis_steps,verbose=False):
    print("checking status of jobs...\n")
    spec_yaml = ConnectionManager.execute_command(ssh_client, 
            "cat /shared/workspace/Pipelines/config/{}/{}_{}.yaml".format(pipeline, pipeline, workflow))
    spec_yaml = yaml.load(spec_yaml)
    
    possible_steps = get_possible_steps(analysis_steps, spec_yaml)

    if verbose:
        print("Your project will go through the following steps:\n\t{}\n".format(
            ", ".join(possible_steps)))

    all_possible_job_dict = get_job_names("all", possible_steps, spec_yaml)
    job_dict = get_job_names(step_name, possible_steps, spec_yaml)

    qstat = ConnectionManager.execute_command(ssh_client, "qstat")
    current_job = get_current_job(ssh_client, qstat)

    if qstat:
        split_qstat = qstat.splitlines()[2:]
    else:
        split_qstat = []

    curr_time = datetime.datetime.utcnow() #for time running in verbose output
    status_conv = {"qw":"queued", "r":"running", "dr":"being deleted", "t":"being transferred"}
    possible_steps = [all_possible_job_dict[x] for x in possible_steps if not x == "done"]
    possible_steps.append("done")
    is_done = True

    for step, job_name in job_dict.items():
        if verbose:
            print("The {} step calls the {} script on the cluster".format(step, job_name))

        qstat_j = get_qstat_j(ssh_client, job_name)

        if "Following jobs do not exist" in qstat_j: #only happens when qstat -j job_name fails
            if not job_name in possible_steps:
                print("The {} step was not specified as a step in analysis_steps".format(step))
                is_done = is_done and True
            elif possible_steps.index(current_job) < possible_steps.index(job_name):
                print("The {} step has not started yet.".format(step))
                is_done = False
            elif check_step_failed(ssh_client, pipeline, workflow, project_name, job_name):
                print("The {} step has finished running, but has failed".format(step))
                print("\tPlease check the logs")
                is_done = is_done and True
            elif check_step_passed(ssh_client, pipeline, workflow, project_name, job_name):
                print("The {} step has finished running without failure".format(step))
                is_done = is_done and True
            else:
                print("The {} step has not started yet.".format(step))
                is_done = False
        else: #job must be in qstat
            print("The {} step is being executed".format(step))

            num_jobs = {x:0 for x in status_conv.keys()}
            job_specifics = defaultdict(list)

            for line in split_qstat:
                line = line.split()
                line_id = line[0]

                #get job name from qstat -j otherwise job name could be longer than qstat allows
                line_job = re.search("job_name:\s+(.*)\n",get_qstat_j(ssh_client, job_name)).group(1)
                if not line_job == job_name:
                    continue

                line_status = line[4]
                num_jobs[line_status] += 1

                month,day,year = map(int,line[5].split("/"))
                hour,minute,second = map(int,line[6].split(":"))
                start_time = datetime.datetime(year,month,day,hour,minute,second)

                line_num_cores = line[-1]
                job_specifics[line_status].append("\tone is currently {} using {} core(s) and was submitted {} ago".format(
                    status_conv[line_status], line_num_cores, _format_timedelta(curr_time-start_time)))

            step_info = {stat:(num,job_specifics[stat]) for stat,num in num_jobs.items()}

            for stat,info_tuple in step_info.items():
                num, details = info_tuple
                if num == 0:
                    continue
                print("There are {} instances of the {} step currently {}".format(num, step, status_conv[stat]))
                if verbose:
                    for det in details:
                        print(det)
        print()

    if current_job == "done" and is_done:
        print("\nYour pipeline has finished")
    print()

def get_possible_steps(analysis_steps, spec_yaml):
    possible_steps = filter(lambda x:x in analysis_steps, spec_yaml["steps"])
    return list(possible_steps)


def get_job_names(step_name, possible_steps, spec_yaml):
    if not step_name == "all":
        return {step_name:_step_to_job(spec_yaml, step_name)}
    else:
        return {step:_step_to_job(spec_yaml, step) for step in possible_steps if not step == "done"}


def get_current_job(ssh_client, qstat):
    if not qstat:
        return "done"
    else:
        qstat_j = get_qstat_j(ssh_client, qstat.splitlines()[2].split()[0])
        for line in qstat_j.splitlines():
            if line.startswith("job_name"):
                return line.split()[1]

def get_qstat_j(ssh_client, job_name):
    return ConnectionManager.execute_command(ssh_client, "qstat -j {}".format(job_name))


def check_step_failed(ssh_client, pipeline, workflow, project_name, job_name):
    ConnectionManager.execute_command(ssh_client, "touch /shared/workspace/logs/{}/{}/{}/*/status.log".format(pipeline,workflow,project_name,job_name))
    status_log_checker = "ls /shared/workspace/logs/{}/{}/{}/*/status.log | xargs grep \"{}.*failed\" | wc -l"

    #reports if step finished and failed or finished and passed
    if int(ConnectionManager.execute_command(ssh_client,
        status_log_checker.format(pipeline,workflow,project_name,job_name))):
        return True

    return False

def check_step_passed(ssh_client, pipeline, workflow, project_name, job_name):
    status_log_checker = "ls /shared/workspace/logs/{}/{}/{}/{}status.log | xargs grep \"{}.*passed\" | wc -l"

    #reports if step finished and failed or finished and passed
    if int(ConnectionManager.execute_command(ssh_client,
        status_log_checker.format(pipeline,workflow,project_name,"*/",job_name))):
        return True
    elif int(ConnectionManager.execute_command(ssh_client,
        status_log_checker.format(pipeline,workflow,project_name,"",job_name)).split('\n')[0]):
        return True

    return False

def _step_to_job(spec_yaml, step_name):
    if step_name in spec_yaml:
        return spec_yaml[step_name]["script_path"].split("/")[-1] + ".sh"
    else:
        raise ValueError("{} is not a possible step in your pipeline.".format(step_name))


def _format_timedelta(td):
    return "{} days, {} hours, and {} minutes".format(td.days, td.seconds//3600, (td.seconds//60)%60)
