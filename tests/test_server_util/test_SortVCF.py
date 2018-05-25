import unittest
import sys
import os
#sys.path.append(os.getcwd().replace("test", "src"))
sys.path.append(os.getcwd().replace("test", "src") + "/cirrus_ngs/server/Pipelines/util")
import SortVCF as svcf
import tempfile 
from getopt import GetoptError


class Tests(unittest.TestCase):

    usage_msg = """
usage:
    python SortVCF.py <vcf_file> <chrom_order> [options]

args:
    vcf_file        path to vcf file to be sorted
                        e.g. /scratch/project/raw.vcf
    chrom_order     space-separated string with target chromosome ordering
                        e.g. "1 2 3" or "chr10 chr2 chr5"
                        the target chromsomes must match chromosome names in vcf file

options:
    -o              specify output file, default stdout
    -h              print this usage message and exit
"""

    dummy_vcf_text = """##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
20     14370   rs6054257 G      A       29   PASS   NS=3;DP=14;AF=0.5;DB;H2           GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
12     17330   .         T      A       3    q10    NS=3;DP=11;AF=0.017               GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3
20     1110696 rs6040355 A      G,T     67   PASS   NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27 2|1:2:0:18,2   2/2:35:4
12     1230237 .         T      .       47   PASS   NS=3;DP=13;AA=T                   GT:GQ:DP:HQ 0|0:54:7:56,60 0|0:48:4:51,51 0/0:61:2
3     1234567 microsat1 GTCT   G,GTACT 50   PASS   NS=3;DP=9;AA=G                    GT:GQ:DP    0/1:35:4       0/2:17:2       1/1:40:3"""

    def test_get_options_help(self):
        args = ["path", "1 2 3", "-o", "hello", "-h"]
        with self.assertRaises(SystemExit) as err:
            svcf.get_options(args)
        self.assertEqual(err.exception.code, self.usage_msg)

    def test_get_options_fake_option(self):
        args = ["-c", "-o", "hello", "-h"]
        with self.assertRaises(GetoptError) as err:
            svcf.get_options(args)

        args = ["-c", "-o", "-h"]
        with self.assertRaises(GetoptError) as err:
            svcf.get_options(args)

    def test_get_options_missing_args(self):
        no_file_args = ["2 3 4", "-o", "outfile"]
        with self.assertRaises(SystemExit) as err:
            svcf.get_options(no_file_args)
        self.assertEqual(err.exception.code, self.usage_msg)

        no_chrom_order_args = ["file", "-o", "outfile"]
        with self.assertRaises(SystemExit) as err:
            svcf.get_options(no_chrom_order_args)
        self.assertEqual(err.exception.code, self.usage_msg)

    def test_get_options_correct_usage(self):
        args = ["path/to/vcf", "1 2 3 4", "-o", "./out"]
        vcf_file, chrom_order, out_file = svcf.get_options(args)

        self.assertEqual(vcf_file, "path/to/vcf")
        self.assertEqual(chrom_order, ["1","2","3","4"])
        self.assertEqual(out_file.name, "./out")

        args = ["path/to/vcf", "1 2 3 4"]
        vcf_file, chrom_order, out_file = svcf.get_options(args)

        self.assertEqual(vcf_file, "path/to/vcf")
        self.assertEqual(chrom_order, ["1","2","3","4"])
        self.assertEqual(out_file, sys.stdout)

    def test_print_headers(self):
        vcf_file = tempfile.NamedTemporaryFile(mode="w+", delete=False)
        vcf_file.write(self.dummy_vcf_text)
        vcf_file.close()

        out_file = tempfile.NamedTemporaryFile(mode="w+", delete=False)

        headers = """##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
"""
        correct_offset = len(headers)

        self.assertEqual(correct_offset, svcf.print_headers(vcf_file.name, out_file))
        out_file.close()

        with open(out_file.name) as f:
            self.assertEqual(headers ,f.read())

    def test_parse_variants(self):
        vcf_file = tempfile.NamedTemporaryFile(mode="w+", delete=False)
        vcf_file.write(self.dummy_vcf_text)
        vcf_file.close()
        out_file = tempfile.NamedTemporaryFile(mode="w+")

        header_offset = svcf.print_headers(vcf_file.name, out_file)
        out_file.close()
        chrom_dict = svcf.parse_variants(vcf_file.name, ["3", "20", "12"], header_offset)

        entry_offsets = [1060, 1201, 1338, 1475, 1612]
        correct_dict = {"3":[(1234567, 1612)], "12":[(17330, 1201), (1230237, 1475)], "20":[(14370,1060), (1110696, 1338)]}
        self.assertEqual(chrom_dict, correct_dict)

        #makes sure that error is raised when chrom not in chromosome order is found
        out_file = tempfile.NamedTemporaryFile(mode="w+")
        with self.assertRaises(SystemExit) as err:
            svcf.parse_variants(vcf_file.name, ["5"], header_offset)
        self.assertEqual(err.exception.code, """
ERROR
    Found chromosome in vcf_file not found in chrom_order input
        chrom_order = \"{}\"
        failing chromosome = \"{}\"
""".format("5", "20") + self.usage_msg)

    def test_sort_vcf(self):
        vcf_file = tempfile.NamedTemporaryFile(mode="w+", delete=False)
        vcf_file.write(self.dummy_vcf_text)
        vcf_file.close()
        out_file = tempfile.NamedTemporaryFile(mode="w+", delete=False)

        header_offset = svcf.print_headers(vcf_file.name, out_file)
        chrom_dict = svcf.parse_variants(vcf_file.name, ["3", "20", "12"], header_offset)

        svcf.sort_vcf(vcf_file.name, chrom_dict, ["3", "20", "12"], out_file)
        out_file.close()

        self.maxDiff = None

        with open(out_file.name) as f:
            self.assertEqual(f.read(), """##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
3     1234567 microsat1 GTCT   G,GTACT 50   PASS   NS=3;DP=9;AA=G                    GT:GQ:DP    0/1:35:4       0/2:17:2       1/1:40:3
20     14370   rs6054257 G      A       29   PASS   NS=3;DP=14;AF=0.5;DB;H2           GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
20     1110696 rs6040355 A      G,T     67   PASS   NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27 2|1:2:0:18,2   2/2:35:4
12     17330   .         T      A       3    q10    NS=3;DP=11;AF=0.017               GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3
12     1230237 .         T      .       47   PASS   NS=3;DP=13;AA=T                   GT:GQ:DP:HQ 0|0:54:7:56,60 0|0:48:4:51,51 0/0:61:2""")















        








if __name__ == "__main__":
    unittest.main()
