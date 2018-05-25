import unittest
import sys
import os
import cirrusngs.util.QstatParser as QstatParser
import datetime

class Tests(unittest.TestCase):
    spec_yaml = {"steps": ["step1", "step2", "step3", "step4"], 
            "step1": { "script_path": "test1", 
                "download_suffix": None, 
                "input_is_output": False, 
                "can_be_zipped": True, 
                "uses_chromosomes": False }, 
            "step2": { "script_path": "hello/test2", 
                "download_suffix": None, 
                "input_is_output": False, 
                "can_be_zipped": True, 
                "uses_chromosomes": False }, 
            "step3": { "script_path": "extra/stuff/here/test3", 
                "download_suffix": None, 
                "input_is_output": False, 
                "can_be_zipped": True, 
                "uses_chromosomes": False }, 
            "step4": { "script_path": "test4", 
                "download_suffix": None, 
                "input_is_output": False, 
                "can_be_zipped": True, 
                "uses_chromosomes": False }}

#    def test_get_current_job(self):
#        qstat = ""
#        output = QstatParser.get_current_job(qstat)
#        self.assertEqual(output, "done")
#        qstat = """job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID
#-----------------------------------------------------------------------------------------------------------------
#362 0.00000 slow.sh    ec2-user     qw    07/07/2017 18:55:45                                    1
#407 0.50000 slow.sh    user         r     07/07/2017 18:23:12     queue_here                     1
#"""
#
#        output = QstatParser.get_current_job(qstat)
#        self.assertEqual(output, "slow.sh")
#
#        qstat = """job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID
#-----------------------------------------------------------------------------------------------------------------
#362 0.00000 new.sh    ec2-user     qw    07/07/2017 18:55:45                                    1
#407 0.50000 new.sh    user         r     07/07/2017 18:23:12     queue_here                     1
#"""
#        output = QstatParser.get_current_job(qstat)
#        self.assertEqual(output, "new.sh")


    def test_format_timedelta(self):
        start = datetime.datetime(2000,10,10,hour=1)
        end = datetime.datetime(2000,10,10,hour=2)
        td = end - start
        output = QstatParser._format_timedelta(td)
        self.assertEqual(output, "0 days, 1 hours, and 0 minutes")

        start = datetime.datetime(2000,10,10,hour=1, minute=20)
        end = datetime.datetime(2000,10,10,hour=2, minute=10)
        td = end - start
        output = QstatParser._format_timedelta(td)
        self.assertEqual(output, "0 days, 0 hours, and 50 minutes")

        start = datetime.datetime(2000,10,9,hour=1, minute=20)
        end = datetime.datetime(2000,10,10,hour=2, minute=10)
        td = end - start
        output = QstatParser._format_timedelta(td)
        self.assertEqual(output, "1 days, 0 hours, and 50 minutes")

    def test_get_possible_steps(self):
        analysis_steps = {"step1", "step3"}
        output = QstatParser.get_possible_steps(analysis_steps, Tests.spec_yaml)
        self.assertEqual(output, ["step1", "step3"])
        
        analysis_steps = {"step1", "trash", "step3"}
        output = QstatParser.get_possible_steps(analysis_steps, Tests.spec_yaml)
        self.assertEqual(output, ["step1", "step3"])

        analysis_steps = {"step1", "step3", "step4", "step2"}
        output = QstatParser.get_possible_steps(analysis_steps, Tests.spec_yaml)
        self.assertEqual(output, ["step1", "step2", "step3", "step4"])

    def test_step_to_job(self):
        steps = ["step1", "step3", "step4"]
        jobs = ["test1.sh", "test3.sh", "test4.sh"]

        for step, job in zip(steps, jobs):
            output = QstatParser._step_to_job(Tests.spec_yaml, step)
            self.assertEqual(output, job)

        with self.assertRaises(ValueError) as err:
            output = QstatParser._step_to_job({"fake":"hello"}, "step1")

    def test_get_job_names(self):
        possible_steps = "step1 step2 step3 step4".split()

        step_name = "all"
        output = QstatParser.get_job_names(step_name, possible_steps, Tests.spec_yaml)
        self.assertEqual(output, {"step1":"test1.sh", "step2":"test2.sh", "step3":"test3.sh", "step4":"test4.sh"})

        steps = ["step1", "step2", "step3", "step4"]
        jobs = ["test1.sh", "test2.sh", "test3.sh", "test4.sh"]
        for step, job in zip(steps, jobs):
            output = QstatParser.get_job_names(step, possible_steps, Tests.spec_yaml)
            self.assertEqual(output, {step:job})


        


    
if __name__ == "__main__":
    unittest.main()
