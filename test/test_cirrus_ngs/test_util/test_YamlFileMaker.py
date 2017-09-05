import unittest
import sys
import os
sys.path.append(os.getcwd().replace("test", "src"))
import cirrus_ngs.util.YamlFileMaker as yfm
import tempfile 
import time

class Tests(unittest.TestCase):

    def test_make_yaml_file_single_end_no_pairs_list(self):
        yaml_file = tempfile.NamedTemporaryFile("w", delete=False)
        yaml_file.close()

        project_name = "test_project"
        analysis_steps = ["fastqc", "trim", "bwa"]
        s3_input_files_address = "s3://fake/s3/path"
        sample_list = [["sample1.fq"], ["sample2.fq"]]
        group_list = ["groupA", "groupB"]
        pair_list = None
        s3_output_files_address = "s3://fake/s3/out/path"
        genome = "hg19"
        style = "histone"

        yfm.make_yaml_file(yaml_file.name, project_name, analysis_steps,
                s3_input_files_address, sample_list, group_list,
                s3_output_files_address, genome, style, None)

        date = time.strftime("%Y/%m/%d")

        correct_output = """project: test_project
analysis:
  - fastqc
  - trim
  - bwa
date: """ + date +"""
upload: s3://fake/s3/out/path
download: s3://fake/s3/path
genome: hg19
style: histone
pairs:
  {}
sample:
- filename: sample1.fq
  description: sample1
  group: groupA
- filename: sample2.fq
  description: sample2
  group: groupB
"""

        temp = open(yaml_file.name, "r")
        yaml_file_contents = temp.read()
        temp.close()

        self.assertEqual(yaml_file_contents, correct_output)

    def test_make_yaml_file_single_end_pairs_list(self):
        yaml_file = tempfile.NamedTemporaryFile("w", delete=False)
        yaml_file.close()

        project_name = "test_project"
        analysis_steps = ["trim", "bwa"]
        s3_input_files_address = "s3://fake/s3/path"
        sample_list = [["sample1.fq"], ["sample1b.fq"], ["sample2.fq"], ["sample2b.fq"]]
        group_list = ["groupA", "groupA", "groupB", "groupB"]
        pair_list = { "sample1": "sample1b", "sample2": "sample2b" }
        s3_output_files_address = "s3://fake/s3/out/path"
        genome = "hg19"
        style = "histone"

        yfm.make_yaml_file(yaml_file.name, project_name, analysis_steps,
                s3_input_files_address, sample_list, group_list,
                s3_output_files_address, genome, style, pair_list)

        date = time.strftime("%Y/%m/%d")

        correct_output = """project: test_project
analysis:
  - trim
  - bwa
date: """ + date +"""
upload: s3://fake/s3/out/path
download: s3://fake/s3/path
genome: hg19
style: histone
pairs:
  sample1: sample1b
  sample2: sample2b
sample:
- filename: sample1.fq
  description: sample1
  group: groupA
- filename: sample1b.fq
  description: sample1b
  group: groupA
- filename: sample2.fq
  description: sample2
  group: groupB
- filename: sample2b.fq
  description: sample2b
  group: groupB
"""
        temp = open(yaml_file.name, "r")
        yaml_file_contents = temp.read()
        temp.close()

        self.assertEqual(yaml_file_contents, correct_output)

    def test_make_yaml_file_paired_end_no_pairs_list(self):
        yaml_file = tempfile.NamedTemporaryFile("w", delete=False)
        yaml_file.close()

        project_name = "test_project"
        analysis_steps = ["bwa"]
        s3_input_files_address = "s3://fake/s3/path"
        sample_list = [["sample1.fq", "reverse1.fq"], ["sample2.fq", "reverse2.fq"]]
        group_list = ["groupA", "groupB"]
        pair_list = None
        s3_output_files_address = "s3://fake/s3/out/path"
        genome = "hg19"
        style = "histone"

        yfm.make_yaml_file(yaml_file.name, project_name, analysis_steps,
                s3_input_files_address, sample_list, group_list,
                s3_output_files_address, genome, style, pair_list)

        date = time.strftime("%Y/%m/%d")

        correct_output = """project: test_project
analysis:
  - bwa
date: """ + date +"""
upload: s3://fake/s3/out/path
download: s3://fake/s3/path
genome: hg19
style: histone
pairs:
  {}
sample:
- filename: sample1.fq, reverse1.fq
  description: sample1
  group: groupA
- filename: sample2.fq, reverse2.fq
  description: sample2
  group: groupB
"""
        temp = open(yaml_file.name, "r")
        yaml_file_contents = temp.read()
        temp.close()

        self.assertEqual(yaml_file_contents, correct_output)

    def test_make_yaml_file_paired_end_pairs_list(self):
        yaml_file = tempfile.NamedTemporaryFile("w", delete=False)
        yaml_file.close()

        project_name = "test_project"
        analysis_steps = ["bwa"]
        s3_input_files_address = "s3://fake/s3/path"
        sample_list = [["sample1.fq", "reverse1.fq"], ["sample2.fq", "reverse2.fq"]]
        group_list = ["groupA", "groupB"]
        pair_list = {"sample1": "sample2"}
        s3_output_files_address = "s3://fake/s3/out/path"
        genome = "hg19"
        style = "histone"

        yfm.make_yaml_file(yaml_file.name, project_name, analysis_steps,
                s3_input_files_address, sample_list, group_list,
                s3_output_files_address, genome, style, pair_list)

        date = time.strftime("%Y/%m/%d")

        correct_output = """project: test_project
analysis:
  - bwa
date: """ + date +"""
upload: s3://fake/s3/out/path
download: s3://fake/s3/path
genome: hg19
style: histone
pairs:
  sample1: sample2
sample:
- filename: sample1.fq, reverse1.fq
  description: sample1
  group: groupA
- filename: sample2.fq, reverse2.fq
  description: sample2
  group: groupB
"""
        temp = open(yaml_file.name, "r")
        yaml_file_contents = temp.read()
        temp.close()

        self.assertEqual(yaml_file_contents, correct_output)

if __name__ == "__main__":
    unittest.main()
