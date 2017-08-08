import unittest
import sys
import os
sys.path.append(os.getcwd().replace("test", "src"))
import cirrus_ngs.util.DesignFileLoader as dfl
import tempfile 

class Tests(unittest.TestCase):
    def test_load_design_file_should_ignore_comments(self):
        design_file = tempfile.NamedTemporaryFile("w+", delete=False)
        design_file.write("##ignore this message\n")
        design_file.write("sample1.fastq\tgroupA\n")
        design_file.write("##ignore this too\n")
        design_file.write("sample2.fq\tgroupA\n")
        design_file.write("##ignore this too\n")
        design_file.write("##ignore this too\n")
        design_file.write("##ignore this too\n")
        design_file.close()

        sample_list, group_list = dfl.load_design_file(design_file.name)

        self.assertEqual(sample_list, [["sample1.fastq"], ["sample2.fq"]])
        self.assertEqual(group_list, ["groupA", "groupA"])

    def test_load_design_file_paired_end_samples(self):
        design_file = tempfile.NamedTemporaryFile("w+", delete=False)
        design_file.write("sample1.fastq,sample1_r1.fq\tgroupA\n")
        design_file.write("sample2.fq,sample2_r2.fastq\tgroupA\n")
        design_file.close()

        sample_list, group_list = dfl.load_design_file(design_file.name)

        self.assertEqual(sample_list, [["sample1.fastq", "sample1_r1.fq"],["sample2.fq", "sample2_r2.fastq"]])
        self.assertEqual(group_list, ["groupA", "groupA"])

    def test_load_design_file_paired_and_single_end_samples(self):
        design_file = tempfile.NamedTemporaryFile("w+", delete=False)
        design_file.write("sample1.fastq,sample1_r1.fq\tgroupA\n")
        design_file.write("sample2.fq\tgroupA\n")
        design_file.close()

        sample_list, group_list = dfl.load_design_file(design_file.name)

        self.assertEqual(sample_list, [["sample1.fastq", "sample1_r1.fq"],["sample2.fq"]])
        self.assertEqual(group_list, ["groupA", "groupA"])

    def test_load_design_file_group_list_order(self):
        design_file = tempfile.NamedTemporaryFile("w+", delete=False)
        design_file.write("sample1.fastq,sample1_r1.fq\tgroupA\n")
        design_file.write("sample2.fq\tgroupB\n")
        design_file.write("sample3.fastq,sample3_r1.fq\tgroupA\n")
        design_file.write("samp5.fastq\tgroupC\n")
        design_file.close()

        sample_list, group_list = dfl.load_design_file(design_file.name)

        self.assertEqual(group_list, ["groupA", "groupB", "groupA", "groupC"])
        self.assertEqual(sample_list, [["sample1.fastq", "sample1_r1.fq"], ["sample2.fq"], ["sample3.fastq", "sample3_r1.fq"], ["samp5.fastq"]])

    def test_load_design_file_rows_arent_tab_separated(self):
        design_file = tempfile.NamedTemporaryFile("w+", delete=False)
        design_file.write("sample1.fastq,sample1_r1.fq groupA\n")
        design_file.close()

        self.assertRaises(IndexError, dfl.load_design_file, design_file.name)

if __name__ == "__main__":
    unittest.main()
