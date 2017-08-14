import unittest
import sys
import os
sys.path.append("../src/cirrus_ngs/server/Pipelines/")
import WGSPipeline as wp

class Tests(unittest.TestCase):
    def test_separate_file_suffix_unzipped_fastq(self):
        sample_file = "Sample1_R1.fastq"
        file_prefix, file_suffix, is_zipped = wp._separate_file_suffix(sample_file)
        self.assertEqual(file_prefix, "Sample1_R1")
        self.assertEqual(file_suffix, ".fastq")
        self.assertEqual(is_zipped, "False")

        sample_file = "Sample1_R1.trim.garbage.bye.fastq"
        file_prefix, file_suffix, is_zipped = wp._separate_file_suffix(sample_file)
        self.assertEqual(file_prefix, "Sample1_R1.trim.garbage.bye")
        self.assertEqual(file_suffix, ".fastq")
        self.assertEqual(is_zipped, "False")

        sample_file = "Sample1_R1.trim.garbage.fastq.bye"
        file_prefix, file_suffix, is_zipped = wp._separate_file_suffix(sample_file)
        self.assertEqual(file_prefix, "Sample1_R1.trim.garbage")
        self.assertEqual(file_suffix, ".fastq.bye")
        self.assertEqual(is_zipped, "False")

    def test_separate_file_suffix_zipped_fastq(self):
        sample_file = "Sample1_R1.fastq.gz"
        file_prefix, file_suffix, is_zipped = wp._separate_file_suffix(sample_file)
        self.assertEqual(file_prefix, "Sample1_R1")
        self.assertEqual(file_suffix, ".fastq")
        self.assertEqual(is_zipped, "True")

        sample_file = "Sample1_R1.trim.garbage.bye.fastq.gz"
        file_prefix, file_suffix, is_zipped = wp._separate_file_suffix(sample_file)
        self.assertEqual(file_prefix, "Sample1_R1.trim.garbage.bye")
        self.assertEqual(file_suffix, ".fastq")
        self.assertEqual(is_zipped, "True")

        sample_file = "Sample1_R1.trim.garbage.fastq.gz.bye"
        file_prefix, file_suffix, is_zipped = wp._separate_file_suffix(sample_file)
        self.assertEqual(file_prefix, "Sample1_R1.trim.garbage")
        self.assertEqual(file_suffix, ".fastq.bye")
        self.assertEqual(is_zipped, "True")

    def test_separate_file_suffix_unzipped_fq(self):
        sample_file = "Sample1_R1.fq"
        file_prefix, file_suffix, is_zipped = wp._separate_file_suffix(sample_file)
        self.assertEqual(file_prefix, "Sample1_R1")
        self.assertEqual(file_suffix, ".fq")
        self.assertEqual(is_zipped, "False")

        sample_file = "Sample1_R1.trim.garbage.bye.fq"
        file_prefix, file_suffix, is_zipped = wp._separate_file_suffix(sample_file)
        self.assertEqual(file_prefix, "Sample1_R1.trim.garbage.bye")
        self.assertEqual(file_suffix, ".fq")
        self.assertEqual(is_zipped, "False")

        sample_file = "Sample1_R1.trim.garbage.fq.bye"
        file_prefix, file_suffix, is_zipped = wp._separate_file_suffix(sample_file)
        self.assertEqual(file_prefix, "Sample1_R1.trim.garbage")
        self.assertEqual(file_suffix, ".fq.bye")
        self.assertEqual(is_zipped, "False")

    def test_separate_file_suffix_zipped_fq(self):
        sample_file = "Sample1_R1.fq.gz"
        file_prefix, file_suffix, is_zipped = wp._separate_file_suffix(sample_file)
        self.assertEqual(file_prefix, "Sample1_R1")
        self.assertEqual(file_suffix, ".fq")
        self.assertEqual(is_zipped, "True")

        sample_file = "Sample1_R1.trim.garbage.bye.fq.gz"
        file_prefix, file_suffix, is_zipped = wp._separate_file_suffix(sample_file)
        self.assertEqual(file_prefix, "Sample1_R1.trim.garbage.bye")
        self.assertEqual(file_suffix, ".fq")
        self.assertEqual(is_zipped, "True")

        sample_file = "Sample1_R1.trim.garbage.fq.gz.bye"
        file_prefix, file_suffix, is_zipped = wp._separate_file_suffix(sample_file)
        self.assertEqual(file_prefix, "Sample1_R1.trim.garbage")
        self.assertEqual(file_suffix, ".fq.bye")
        self.assertEqual(is_zipped, "True")

    def test_sample_argument_generator(self):
        sample_list = [{"filename":"sample1.fq.gz", "download":"s3://fake/s3/dl", "description":"sample1", "group":"groupA"}, 
                {"filename":"sample2.fastq, sample2_rev.fastq", "download":"s3://fake/s3/dl", "description":"sample2", "group":"groupB"}]
        output_address = "s3://fake/address"
        config_dictionary = {"script_name": "fastqc", "download_suffix": None,
                "input_is_output": False, "can_be_zipped": True, "uses_chromosomes":False}

        arg_matrix = []
        for arg_list in wp._sample_argument_generator(sample_list, "s3://fake/output", config_dictionary):
            arg_matrix.append(arg_list)

        self.assertEqual(len(arg_matrix), 2)
        
        for arg_list in arg_matrix:
            self.assertEqual(len(arg_list), 8)

        correct_output = [[".fq", "/scratch", "sample1", "NULL", 
            "s3://fake/s3/dl", "s3://fake/output", 
            "/shared/workspace/logs/DNASeq/", "True"],
            [".fastq", "/scratch", "sample2", "sample2_rev",
                "s3://fake/s3/dl", "s3://fake/output",
                "/shared/workspace/logs/DNASeq/", "False"]]

        self.assertEqual(correct_output, arg_matrix)

        config_dictionary = {"script_name": "fastqc", "download_suffix": None,
                "input_is_output": True, "can_be_zipped": False, "uses_chromosomes":False}

        arg_matrix = []
        for arg_list in wp._sample_argument_generator(sample_list, "s3://fake/output", config_dictionary):
            arg_matrix.append(arg_list)

        self.assertEqual(len(arg_matrix), 2)
        
        for arg_list in arg_matrix:
            self.assertEqual(len(arg_list), 8)

        correct_output = [[".fq", "/scratch", "sample1", "NULL", 
            "s3://fake/output", "s3://fake/output", 
            "/shared/workspace/logs/DNASeq/", "False"],
            [".fastq", "/scratch", "sample2", "sample2_rev",
                "s3://fake/output", "s3://fake/output",
                "/shared/workspace/logs/DNASeq/", "False"]]

        self.assertEqual(correct_output, arg_matrix)

#    def test_by_group_argument_generator(self):
#        group_list = {}
#        output_address ="s3://fake/output"
#        config_dictionary = {}
#        self.fail("this test isn't done yet")




if __name__ == "__main__":
    unittest.main()
