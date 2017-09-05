import unittest
import sys
import os
sys.path.append("../src/cirrus_ngs/server/Pipelines/")
import Pipeline

class Tests(unittest.TestCase):
    # Next 4 tests test the different combinations of fastq/fq and zipped/unzipped extensions
    # _separate_file_suffix should be able to correctly handle any of these combinations
    def test_separate_file_suffix_unzipped_fastq(self):
        sample_file = "Sample1_R1.fastq"
        file_prefix, file_suffix, is_zipped = Pipeline._separate_file_suffix(sample_file)
        self.assertEqual(file_prefix, "Sample1_R1")
        self.assertEqual(file_suffix, ".fastq")
        self.assertEqual(is_zipped, "False")

        sample_file = "Sample1_R1.trim.garbage.bye.fastq"
        file_prefix, file_suffix, is_zipped = Pipeline._separate_file_suffix(sample_file)
        self.assertEqual(file_prefix, "Sample1_R1.trim.garbage.bye")
        self.assertEqual(file_suffix, ".fastq")
        self.assertEqual(is_zipped, "False")

        sample_file = "Sample1_R1.trim.garbage.fastq.bye"
        file_prefix, file_suffix, is_zipped = Pipeline._separate_file_suffix(sample_file)
        self.assertEqual(file_prefix, "Sample1_R1.trim.garbage")
        self.assertEqual(file_suffix, ".fastq.bye")
        self.assertEqual(is_zipped, "False")

    def test_separate_file_suffix_zipped_fastq(self):
        sample_file = "Sample1_R1.fastq.gz"
        file_prefix, file_suffix, is_zipped = Pipeline._separate_file_suffix(sample_file)
        self.assertEqual(file_prefix, "Sample1_R1")
        self.assertEqual(file_suffix, ".fastq")
        self.assertEqual(is_zipped, "True")

        sample_file = "Sample1_R1.trim.garbage.bye.fastq.gz"
        file_prefix, file_suffix, is_zipped = Pipeline._separate_file_suffix(sample_file)
        self.assertEqual(file_prefix, "Sample1_R1.trim.garbage.bye")
        self.assertEqual(file_suffix, ".fastq")
        self.assertEqual(is_zipped, "True")

        sample_file = "Sample1_R1.trim.garbage.fastq.gz.bye"
        file_prefix, file_suffix, is_zipped = Pipeline._separate_file_suffix(sample_file)
        self.assertEqual(file_prefix, "Sample1_R1.trim.garbage")
        self.assertEqual(file_suffix, ".fastq.bye")
        self.assertEqual(is_zipped, "True")

    def test_separate_file_suffix_unzipped_fq(self):
        sample_file = "Sample1_R1.fq"
        file_prefix, file_suffix, is_zipped = Pipeline._separate_file_suffix(sample_file)
        self.assertEqual(file_prefix, "Sample1_R1")
        self.assertEqual(file_suffix, ".fq")
        self.assertEqual(is_zipped, "False")

        sample_file = "Sample1_R1.trim.garbage.bye.fq"
        file_prefix, file_suffix, is_zipped = Pipeline._separate_file_suffix(sample_file)
        self.assertEqual(file_prefix, "Sample1_R1.trim.garbage.bye")
        self.assertEqual(file_suffix, ".fq")
        self.assertEqual(is_zipped, "False")

        sample_file = "Sample1_R1.trim.garbage.fq.bye"
        file_prefix, file_suffix, is_zipped = Pipeline._separate_file_suffix(sample_file)
        self.assertEqual(file_prefix, "Sample1_R1.trim.garbage")
        self.assertEqual(file_suffix, ".fq.bye")
        self.assertEqual(is_zipped, "False")

    def test_separate_file_suffix_zipped_fq(self):
        sample_file = "Sample1_R1.fq.gz"
        file_prefix, file_suffix, is_zipped = Pipeline._separate_file_suffix(sample_file)
        self.assertEqual(file_prefix, "Sample1_R1")
        self.assertEqual(file_suffix, ".fq")
        self.assertEqual(is_zipped, "True")

        sample_file = "Sample1_R1.trim.garbage.bye.fq.gz"
        file_prefix, file_suffix, is_zipped = Pipeline._separate_file_suffix(sample_file)
        self.assertEqual(file_prefix, "Sample1_R1.trim.garbage.bye")
        self.assertEqual(file_suffix, ".fq")
        self.assertEqual(is_zipped, "True")

        sample_file = "Sample1_R1.trim.garbage.fq.gz.bye"
        file_prefix, file_suffix, is_zipped = Pipeline._separate_file_suffix(sample_file)
        self.assertEqual(file_prefix, "Sample1_R1.trim.garbage")
        self.assertEqual(file_suffix, ".fq.bye")
        self.assertEqual(is_zipped, "True")

    def test_by_all_samples_argument_generator(self):
        sample_list = [ {   "filename": "sample1_forward-fq,sample1_backward.fq",
                            "description": "sample1_forward",
                            "group": "groupA"   },
                        {   "filename": "sample2.fq",
                            "description": "sample2",
                            "group": "groupB"   }   ]

        config_dict = { "script_name": "test",
                        "download_suffix": None,
                        "input_is_output": False,
                        "can_be_zipped": True,
                        "uses_chromosomes": False   }
        self.fail()

    def test_sample_argument_generator(self):
        sample_list = [ {   "filename": "sample1_forward.fq,sample1_backward.fq",
                            "description": "sample1_forward",
                            "group": "groupA"   },
                        {   "filename": "sample2.fq.gz",
                            "description": "sample2",
                            "group": "groupB"   }   ]

        config_dicts = [{"script_name": "test",
                        "download_suffix": None,
                        "input_is_output": False,
                        "can_be_zipped": True,
                        "uses_chromosomes": False},
                        {"script_name": "test",
                        "download_suffix": ".ext{}.ext2",
                        "input_is_output": False,
                        "can_be_zipped": False,
                        "uses_chromosomes": False},
                        {"script_name": "test",
                        "download_suffix": ".ext{}.ext2",
                        "input_is_output": True,
                        "can_be_zipped": False,
                        "uses_chromosomes": True}]

        base_correct_result = "test_proj {} /scratch {} {} {} {} /path/to/logs {}"
        curr_correct_result = base_correct_result


        result1 = [curr_correct_result.format(".fq", "sample1_forward", "sample1_backward", "s3://path/to/input", "s3://path/to/output/test_proj/sample1_forward", "False").split() ]
        curr_correct_result = base_correct_result
        result1.append(curr_correct_result.format(".fq", "sample2", "NULL", "s3://path/to/input", "s3://path/to/output/test_proj/sample2", "True").split() )
        curr_correct_result = base_correct_result

        result2 = [curr_correct_result.format(".ext.fq.ext2", "sample1_forward", "sample1_backward", "s3://path/to/input", "s3://path/to/output/test_proj/sample1_forward", "False").split() ] 
        curr_correct_result = base_correct_result
        result2.append(curr_correct_result.format(".ext.fq.ext2", "sample2", "NULL", "s3://path/to/input", "s3://path/to/output/test_proj/sample2", "False").split() )
        curr_correct_result = base_correct_result

        result3 = [curr_correct_result.format(".ext{}.ext2", "sample1_forward", "sample1_backward", "s3://path/to/output/test_proj/sample1_forward",  "s3://path/to/output/test_proj/sample1_forward", "False").split() ] 
        curr_correct_result = base_correct_result
        result3.append(curr_correct_result.format(".ext{}.ext2", "sample2", "NULL", "s3://path/to/output/test_proj/sample2",  "s3://path/to/output/test_proj/sample2", "False").split() )


        correct_results = [result1, result2, result3]
        for index, config in enumerate(config_dicts):
            result = []
            for output in Pipeline._sample_argument_generator("test_proj",
                    sample_list, "s3://path/to/input", "s3://path/to/output",
                    config, "/path/to/logs"):
                result.append(output)

            self.assertEqual(result, correct_results[index])

    def test_by_group_argument_generator(self):
        group_list = {  "groupA": [ ["sample1_R1", ".fq", "False"], ["sample2", ".fq", "False"]],
                        "groupB": [ ["sample3_R1", ".fastq", "True"]] }

        config_dicts = [{"script_name": "test",
                        "download_suffix": None,
                        "input_is_output": False,
                        "can_be_zipped": True,
                        "uses_chromosomes": False},
                        {"script_name": "test",
                        "download_suffix": ".ext{}.ext2",
                        "input_is_output": False,
                        "can_be_zipped": False,
                        "uses_chromosomes": False},
                        {"script_name": "test",
                        "download_suffix": ".ext{}.ext2",
                        "input_is_output": True,
                        "can_be_zipped": False,
                        "uses_chromosomes": True}]

        base_correct_result = "test_proj {} /scratch {} {} {} {} /path/to/logs {}"
        curr_correct_result = base_correct_result


        result1 = [curr_correct_result.format(".fq", "groupA", "NULL", "s3://path/to/input", "s3://path/to/output/test_proj/groupA", "False").split() + ["sample1_R1 sample2"]]
        curr_correct_result = base_correct_result
        result1.append(curr_correct_result.format(".fastq", "groupB", "NULL", "s3://path/to/input", "s3://path/to/output/test_proj/groupB", "True").split() + ["sample3_R1"] )
        curr_correct_result = base_correct_result

        result2 = [curr_correct_result.format(".ext.fq.ext2", "groupA", "NULL", "s3://path/to/input", "s3://path/to/output/test_proj/groupA", "False").split() + ["sample1_R1 sample2"]] 
        curr_correct_result = base_correct_result
        result2.append(curr_correct_result.format(".ext.fastq.ext2", "groupB", "NULL", "s3://path/to/input", "s3://path/to/output/test_proj/groupB", "False").split() + ["sample3_R1"])
        curr_correct_result = base_correct_result

        result3 = [curr_correct_result.format(".ext{}.ext2", "groupA", "NULL", "s3://path/to/output",  "s3://path/to/output/test_proj/groupA", "False").split() + ["sample1_R1 sample2"] ]
        curr_correct_result = base_correct_result
        result3.append(curr_correct_result.format(".ext{}.ext2", "groupB", "NULL", "s3://path/to/output",  "s3://path/to/output/test_proj/groupB", "False").split() + ["sample3_R1"])

        correct_results = [result1, result2, result3]
        for index, config in enumerate(config_dicts):
            result = []
            for output in Pipeline._by_group_argument_generator("test_proj",
                    group_list, "s3://path/to/input", "s3://path/to/output",
                    config, "/path/to/logs"):
                result.append(output)

            self.assertEqual(result, correct_results[index])

    def test_by_pair_argument_generator(self):
        group_list = {  "groupA": [ ["sample1_R1", ".fq", "False"], ["sample2", ".fq", "False"]],
                        "groupB": [ ["sample3_R1", ".fastq", "True"], ["sample4", ".fastq", "True"]] }
        pair_list =  {"sample1_R1": "sample2", "sample3_R1": "sample4"}

        config_dicts = [{"script_name": "test",
                        "download_suffix": None,
                        "input_is_output": False,
                        "can_be_zipped": True,
                        "uses_chromosomes": False},
                        {"script_name": "test",
                        "download_suffix": ".ext{}.ext2",
                        "input_is_output": False,
                        "can_be_zipped": False,
                        "uses_chromosomes": False},
                        {"script_name": "test",
                        "download_suffix": ".ext{}.ext2",
                        "input_is_output": True,
                        "can_be_zipped": False,
                        "uses_chromosomes": True}]

        base_correct_result = "test_proj {} /scratch {} {} {} {} /path/to/logs {}"
        curr_correct_result = base_correct_result

        result1 = [curr_correct_result.format(".fq", "sample1_R1", "sample2", "s3://path/to/input", "s3://path/to/output/test_proj/sample1_R1", "False").split() ]
        curr_correct_result = base_correct_result
        result1.append(curr_correct_result.format(".fastq", "sample3_R1", "sample4", "s3://path/to/input", "s3://path/to/output/test_proj/sample3_R1", "True").split() )
        curr_correct_result = base_correct_result

        result2 = [curr_correct_result.format(".ext.fq.ext2", "sample1_R1", "sample2", "s3://path/to/input", "s3://path/to/output/test_proj/sample1_R1", "False").split() ] 
        curr_correct_result = base_correct_result
        result2.append(curr_correct_result.format(".ext.fastq.ext2", "sample3_R1", "sample4", "s3://path/to/input", "s3://path/to/output/test_proj/sample3_R1", "False").split() )
        curr_correct_result = base_correct_result

        result3 = [curr_correct_result.format(".ext{}.ext2", "sample1_R1", "sample2", "s3://path/to/output",  "s3://path/to/output/test_proj/sample1_R1", "False").split() ] 
        curr_correct_result = base_correct_result
        result3.append(curr_correct_result.format(".ext{}.ext2", "sample3_R1", "sample4", "s3://path/to/output",  "s3://path/to/output/test_proj/sample3_R1", "False").split() )

        correct_results = [result1, result2, result3]
        for index, config in enumerate(config_dicts):
            result = []
            for output in Pipeline._by_pair_argument_generator("test_proj",
                    group_list, pair_list, "s3://path/to/input", "s3://path/to/output",
                    config, "/path/to/logs"):
                result.append(output)

            self.assertEqual(result, correct_results[index])


        
        





if __name__ == "__main__":
    unittest.main()
