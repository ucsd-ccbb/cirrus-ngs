import unittest
import sys
import os
sys.path.append(os.getcwd().replace("test", "src"))
sys.path.append(sys.path[-1] + "/cirrus_ngs")
import cirrus_ngs.dnaSeq.WGSPipelineManager as wpm
import tempfile
import re

class Tests(unittest.TestCase):
    def test_sample_argument_generator(self):
        self.fail()


if __name__ == "__main__":
    unittest.main(module=__name__, buffer=True, exit=False)

