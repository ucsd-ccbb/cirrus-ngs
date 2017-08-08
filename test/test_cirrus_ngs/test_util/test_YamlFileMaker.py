import unittest
import sys
import os
sys.path.append(os.getcwd().replace("test", "src"))
import cirrus_ngs.util.YamlFileMaker as yfm
import tempfile 
import time

class Tests(unittest.TestCase):
    def test_make_yaml_file_proper_number_of_lines(self):
        yaml_file = tempfile.NamedTemporaryFile("w", delete=False)
        yaml_file.close()

        project_name = "test_project"
        analysis_steps = ["fastqc", "trim", "bwa"]
        s3_input_files_address = "s3://fake/s3/path"
        sample_list = [["sample1.fq"], ["sample2.fq"]]
        group_list = ["groupA", "groupB"]
        s3_output_files_address = "s3://fake/s3/out/path"
        genome = "NA"
        style = "NA"

        yfm.make_yaml_file(yaml_file.name, project_name, analysis_steps,
                s3_input_files_address, sample_list, group_list,
                s3_output_files_address, genome, style)

        temp = open(yaml_file.name, "r")
        yaml_file_contents = temp.read().split("\n")
        temp.close()

        self.assertEqual(len(yaml_file_contents), 18)

    def test_make_yaml_file_single_end_sample_output(self):
        yaml_file = tempfile.NamedTemporaryFile("w", delete=False)
        yaml_file.close()

        project_name = "test_project"
        analysis_steps = ["fastqc", "trim", "bwa"]
        s3_input_files_address = "s3://fake/s3/path"
        sample_list = [["sample1.fq"], ["sample2.fq"]]
        group_list = ["groupA", "groupB"]
        s3_output_files_address = "s3://fake/s3/out/path"
        genome = "NA"
        style = "NA"

        yfm.make_yaml_file(yaml_file.name, project_name, analysis_steps,
                s3_input_files_address, sample_list, group_list,
                s3_output_files_address, genome, style)

        temp = open(yaml_file.name, "r")
        yaml_file_contents = temp.read()
        temp.close()

        date = time.strftime("%Y/%m/%d")
        correct_output = """project: test_project
analysis: fastqc, trim, bwa
date: {}
upload: s3://fake/s3/out/path
genome: NA
style: NA
---
sample:
  filename: sample1.fq
  download: s3://fake/s3/path
  description: sample1
  group: groupA
---
sample:
  filename: sample2.fq
  download: s3://fake/s3/path
  description: sample2
  group: groupB""".format(date)

        self.assertEqual(correct_output, yaml_file_contents)

    def test_make_yaml_file_paired_end_sample_output(self):
        yaml_file = tempfile.NamedTemporaryFile("w", delete=False)
        yaml_file.close()

        project_name = "test_project"
        analysis_steps = ["fastqc", "trim", "bwa"]
        s3_input_files_address = "s3://fake/s3/path"
        sample_list = [["sample1.fq", "reverse1.fq"], ["sample2.fq", "reverse2.fq"]]
        group_list = ["groupA", "groupB"]
        s3_output_files_address = "s3://fake/s3/out/path"
        genome = "NA"
        style = "NA"

        yfm.make_yaml_file(yaml_file.name, project_name, analysis_steps,
                s3_input_files_address, sample_list, group_list,
                s3_output_files_address, genome, style)

        temp = open(yaml_file.name, "r")
        yaml_file_contents = temp.read()
        temp.close()

        date = time.strftime("%Y/%m/%d")
        correct_output = """project: test_project
analysis: fastqc, trim, bwa
date: {}
upload: s3://fake/s3/out/path
genome: NA
style: NA
---
sample:
  filename: sample1.fq, reverse1.fq
  download: s3://fake/s3/path
  description: sample1
  group: groupA
---
sample:
  filename: sample2.fq, reverse2.fq
  download: s3://fake/s3/path
  description: sample2
  group: groupB""".format(date)

        self.assertEqual(correct_output, yaml_file_contents)

    def test_make_yaml_file_proper_indentation(self):
        yaml_file = tempfile.NamedTemporaryFile("w", delete=False)
        yaml_file.close()

        project_name = "test_project"
        analysis_steps = ["fastqc", "trim", "bwa"]
        s3_input_files_address = "s3://fake/s3/path"
        sample_list = [["sample1.fq"], ["sample2.fq"]]
        group_list = ["groupA", "groupB"]
        s3_output_files_address = "s3://fake/s3/out/path"
        genome = "NA"
        style = "NA"

        yfm.make_yaml_file(yaml_file.name, project_name, analysis_steps,
                s3_input_files_address, sample_list, group_list,
                s3_output_files_address, genome, style)

        temp = open(yaml_file.name, "r")
        indented_lines = temp.read().split("---\nsample:\n")[1:]
        temp.close()

        for line in indented_lines:
            self.assertEqual(line[:2], "  ")


        sample_list = [["sample1.fq", "reverse"], ["sample2.fq", "reverse2"], ["sample3.fq"], ["sampe4.fq"]]
        group_list = ["groupA", "groupB", "groupA", "groupB"]

        yfm.make_yaml_file(yaml_file.name, project_name, analysis_steps,
                s3_input_files_address, sample_list, group_list,
                s3_output_files_address, genome, style)

        temp = open(yaml_file.name, "r")
        indented_lines = temp.read().split("---\nsample:\n")[1:]
        temp.close()

        for line in indented_lines:
            self.assertEqual(line[:2], "  ")

if __name__ == "__main__":
    unittest.main()
