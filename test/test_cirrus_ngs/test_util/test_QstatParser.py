import unittest
import sys
import os
sys.path.append(os.getcwd().replace("test", "src"))
import cirrus_ngs.util.QstatParser as test_class

class Tests(unittest.TestCase):
    def test_get_job_ids(self):
        qstat = """job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID
-----------------------------------------------------------------------------------------------------------------
362 0.00000 slow.sh    ec2-user     qw    07/07/2017 18:55:45                                    1
407 0.50000 test.sh    user         r     07/07/2017 18:23:12     queue_here                     1"""
        output = test_class.get_job_ids(qstat)
        #tests basic functionality of get_job_ids
        self.assertEqual(output, [[362,"qw"], [407, "r"]])

        bad_qstat = """job-ID  prior  name       notright     incorrect
-------------------------------------------------------------------
hello"""
        #makes sure that the id is an int and the function fails with strings
        self.assertRaises(ValueError, test_class.get_job_ids, bad_qstat)

        bad_qstat = """job-ID  prior  name       notright     incorrect
-------------------------------------------------------------------
1 2
3 ls -l h"""

        #makes sure the function fails with an incorrect qstat input
        self.assertRaises(IndexError, test_class.get_job_ids, bad_qstat)

if __name__ == "__main__":
    unittest.main()
