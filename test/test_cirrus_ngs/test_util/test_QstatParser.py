import unittest
import sys
import os
sys.path.append(os.getcwd().replace("test", "src"))
sys.path.append(os.getcwd().replace("test", "src") + "/cirrus_ngs")
print(sys.path[-1])
import cirrus_ngs.util.QstatParser as QstatParser
import datetime

class Tests(unittest.TestCase):
    def test_get_current_job(self):
        qstat = ""
        output = QstatParser.get_current_job(qstat)
        self.assertEqual(output, "done")
        qstat = """job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID
-----------------------------------------------------------------------------------------------------------------
362 0.00000 slow.sh    ec2-user     qw    07/07/2017 18:55:45                                    1
407 0.50000 slow.sh    user         r     07/07/2017 18:23:12     queue_here                     1
"""

        output = QstatParser.get_current_job(qstat)
        self.assertEqual(output, "slow.sh")

        qstat = """job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID
-----------------------------------------------------------------------------------------------------------------
362 0.00000 new.sh    ec2-user     qw    07/07/2017 18:55:45                                    1
407 0.50000 new.sh    user         r     07/07/2017 18:23:12     queue_here                     1
"""
        output = QstatParser.get_current_job(qstat)
        self.assertEqual(output, "new.sh")


    def test_step_to_job(self):
        config_dicts = [{"script_path": "test",
                        "download_suffix": None,
                        "input_is_output": False,
                        "can_be_zipped": True,
                        "uses_chromosomes": False},
                        {"script_path": "test",
                        "download_suffix": ".ext{}.ext2",
                        "input_is_output": False,
                        "can_be_zipped": False,
                        "uses_chromosomes": False},
                        {"script_path": "test",
                        "download_suffix": ".ext{}.ext2",
                        "input_is_output": True,
                        "can_be_zipped": False,
                        "uses_chromosomes": True}]
        step_name = "test_step"

        for dic in config_dicts:
            output = QstatParser._step_to_job({"test_step":dic}, step_name)
            self.assertEqual(output, "test.sh")
            with self.assertRaises(ValueError) as err:
                output = QstatParser._step_to_job({"fake":dic}, step_name)


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

    
    
    
    
    
    
if __name__ == "__main__":
    unittest.main()
