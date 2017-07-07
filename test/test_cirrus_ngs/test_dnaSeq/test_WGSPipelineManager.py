import unittest
import sys
import os
sys.path.append(os.getcwd().replace("test", "src"))
sys.path.append(sys.path[-1] + "/cirrus_ngs")
import cirrus_ngs.dnaSeq.WGSPipelineManager as WGSPipelineManager
import cirrus_ngs.cfnCluster.ConnectionManager as ConnectionManager

class Tests(unittest.TestCase):
    def test_check_status(self):
        assertTrue(False)


if __name__ == "__main__":
    unittest.main()

