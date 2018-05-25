import unittest
import sys
import os
import cirrusngs.util.DesignFileLoader as dfl
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

        sample_list, group_list, pair_list = dfl.load_design_file(design_file.name)

        self.assertEqual(sample_list, [["sample1.fastq"], ["sample2.fq"]])
        self.assertEqual(group_list, ["groupA", "groupA"])
        self.assertEqual(pair_list, {})

    def test_load_design_file_paired_end_samples_no_pair_list(self):
        design_file = tempfile.NamedTemporaryFile("w+", delete=False)
        design_file.write("sample1.fastq,sample1_r1.fq\tgroupA\n")
        design_file.write("sample2.fq,sample2_r2.fastq\tgroupA\n")
        design_file.close()

        sample_list, group_list, pair_list = dfl.load_design_file(design_file.name)

        self.assertEqual(sample_list, [["sample1.fastq", "sample1_r1.fq"],["sample2.fq", "sample2_r2.fastq"]])
        self.assertEqual(group_list, ["groupA", "groupA"])
        self.assertEqual(pair_list, {})

    def test_load_design_file_paired_end_samples_pair_list(self):
        design_file = tempfile.NamedTemporaryFile("w+", delete=False)
        design_file.write("sample1.fastq,sample1_r1.fq\tgroupA\tNormal\n")
        design_file.write("sample2.fq,sample2_r2.fastq\tgroupA\tTumor\n")
        design_file.close()

        sample_list, group_list, pair_list = dfl.load_design_file(design_file.name)

        self.assertEqual(sample_list, [["sample1.fastq", "sample1_r1.fq"],["sample2.fq", "sample2_r2.fastq"]])
        self.assertEqual(group_list, ["groupA", "groupA"])
        self.assertEqual(pair_list, {"sample1": "sample2"})

    def test_load_design_file_paired_and_single_end_samples_no_pair_list(self):
        design_file = tempfile.NamedTemporaryFile("w+", delete=False)
        design_file.write("sample1.fastq,sample1_r1.fq\tgroupA\n")
        design_file.write("sample2.fq\tgroupA\n")
        design_file.close()

        sample_list, group_list, pair_list = dfl.load_design_file(design_file.name)

        self.assertEqual(sample_list, [["sample1.fastq", "sample1_r1.fq"],["sample2.fq"]])
        self.assertEqual(group_list, ["groupA", "groupA"])
        self.assertEqual(pair_list, {})

    def test_load_design_file_paired_and_single_end_samples_pair_list(self):
        design_file = tempfile.NamedTemporaryFile("w+", delete=False)
        design_file.write("sample1.fastq,sample1_r1.fq\tgroupA\tChip\n")
        design_file.write("sample2.fq\tgroupA\tInput\n")
        design_file.close()

        sample_list, group_list, pair_list = dfl.load_design_file(design_file.name)

        self.assertEqual(sample_list, [["sample1.fastq", "sample1_r1.fq"],["sample2.fq"]])
        self.assertEqual(group_list, ["groupA", "groupA"])
        self.assertEqual(pair_list, {"sample1": "sample2"})

    def test_load_design_file_group_list_order(self):
        design_file = tempfile.NamedTemporaryFile("w+", delete=False)
        design_file.write("sample1.fastq,sample1_r1.fq\tgroupA\n")
        design_file.write("sample2.fq\tgroupB\n")
        design_file.write("sample3.fastq,sample3_r1.fq\tgroupA\n")
        design_file.write("samp5.fastq\tgroupC\n")
        design_file.close()

        sample_list, group_list, pair_list = dfl.load_design_file(design_file.name)

        self.assertEqual(group_list, ["groupA", "groupB", "groupA", "groupC"])
        self.assertEqual(sample_list, [["sample1.fastq", "sample1_r1.fq"], ["sample2.fq"], ["sample3.fastq", "sample3_r1.fq"], ["samp5.fastq"]])
        self.assertEqual(pair_list, {})

    def test_load_design_file_pairs_out_of_order(self):
        design_file = tempfile.NamedTemporaryFile("w+", delete=False)
        design_file.write("sample2.fq\tgroupA\tInput\n")
        design_file.write("sample1.fastq,sample1_r1.fq\tgroupA\tChip\n")
        design_file.close()

        sample_list, group_list, pair_list = dfl.load_design_file(design_file.name)

        self.assertEqual(sample_list, [["sample2.fq"], ["sample1.fastq", "sample1_r1.fq"]])
        self.assertEqual(group_list, ["groupA", "groupA"])
        self.assertEqual(pair_list, {"sample1": "sample2"})

    def test_load_design_file_pairs_not_next_to_each_other(self):
        design_file = tempfile.NamedTemporaryFile("w+", delete=False)
        design_file.write("sample1.fastq,sample1_r1.fq\tgroupA\tChip\n")
        design_file.write("sample3.fq\tgroupB\tChip\n")
        design_file.write("sample2.fq\tgroupA\tInput\n")
        design_file.write("sample4.fq\tgroupB\tInput\n")
        design_file.close()

        sample_list, group_list, pair_list = dfl.load_design_file(design_file.name)

        self.assertEqual(sample_list, [["sample1.fastq", "sample1_r1.fq"], ["sample3.fq"], ["sample2.fq"], ["sample4.fq"]])
        self.assertEqual(group_list, ["groupA", "groupB", "groupA", "groupB"])
        self.assertEqual(pair_list, {"sample1": "sample2", "sample3": "sample4"})


    def test_load_design_file_rows_arent_tab_separated(self):
        design_file = tempfile.NamedTemporaryFile("w+", delete=False)
        design_file.write("sample1.fastq,sample1_r1.fq groupA\n")
        design_file.close()

        self.assertRaises(IndexError, dfl.load_design_file, design_file.name)

    def test_load_design_file_third_field_has_wrong_value(self):
        #third field should be Normal, Tumor, Chip, or Input
        design_file = tempfile.NamedTemporaryFile("w+", delete=False)
        design_file.write("sample1.fastq,sample1_r1.fq\tgroupA\tHello\n")
        design_file.close()

        self.assertRaises(ValueError, dfl.load_design_file, design_file.name)

if __name__ == "__main__":
    unittest.main()
