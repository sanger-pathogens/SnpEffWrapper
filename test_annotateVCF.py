#!/usr/bin/env python3

from io import StringIO
from unittest.mock import patch, MagicMock

from annotateVCF import *

class TestAnnotateVCF(unittest.TestCase):
  def setUp(self):
    class FakeArgs(object):
      pass
    self.fake_args = FakeArgs()

  @patch('annotateVCF.argparse.ArgumentParser')
  @patch('annotateVCF.shutil')
  def test_snpeff_not_in_path(self, shutil_mock, argument_parser_mock):
    parsed_args = MagicMock()
    parsed_args.snpeff_exec = 'foobar'
    argument_parser = MagicMock()
    argument_parser.parse_args.return_value = parsed_args
    argument_parser_mock.return_value = argument_parser

    actual_args = parse_arguments()
    self.assertEqual(actual_args.snpeff_exec, 'foobar')

    parsed_args.snpeff_exec = None
    shutil_mock.which.return_value = None
    self.assertRaises(MissingSNPEffError, parse_arguments)
    shutil_mock.which.assert_any_call('snpeff')
    shutil_mock.which.reset_mock()

    parsed_args.snpeff_exec = None
    shutil_mock.which.return_value = '/bin/snpEff'
    actual_args = parse_arguments()
    shutil_mock.which.assert_any_call('snpeff')
    self.assertEqual(actual_args.snpeff_exec, '/bin/snpEff')

  def test_parse_coding_table(self):
    coding_table_str = 'default: Bacterial_and_Plant_Plastid'
    expected = {'default': 'Bacterial_and_Plant_Plastid'}
    actual = parse_coding_table(coding_table_str)
    self.assertEqual(actual, expected)

    coding_table_str = 'foo: bar'
    expected = {'foo': 'bar'}
    actual = parse_coding_table(coding_table_str)
    self.assertEqual(actual, expected)

    coding_table_str = '{foo: bar, default: Bacterial_and_Plant_Plastid}'
    expected = {'foo': 'bar', 'default': 'Bacterial_and_Plant_Plastid'}
    actual = parse_coding_table(coding_table_str)
    self.assertEqual(actual, expected)

  def test_get_gff_contigs(self):
    fake_gff = StringIO("""\
##gff-version 3
##sequence-region CHROM1 1 4000000
##sequence-region PLASMID1 1 2000000
##sequence-region PLASMID2 1 1000000
CHROM1	EMBL	databank_entry	1	4000000	.	+	.	Some=Data
CHROM1	EMBL	CDS	100	400	.	+	0	Some=More	data
PLASMID1	EMBL	databank_entry	1	2000000	.	+	.	Some=Plasmid	data
PLASMID1	EMBL	CDS	300	700	.	+	0	Some=More	plasmid	data
""")
    expected_contigs = ['CHROM1', 'PLASMID1']
    actual_contigs = get_gff_contigs(fake_gff)
    self.assertEqual(actual_contigs, expected_contigs)

  def test_get_vcf_contigs(self):
    fake_vcf = StringIO("""\
##fileformat=VCFv4.1
##contig=<ID=1,length=4000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	10000_2#3	10000_2#4
CHROM1	100	.	G	A	.	.	.	GT	0	1
CHROM1	200	.	G	A	.	.	.	GT	0	1
PLASMID1	50	.	G	T	.	.	.	GT	1	0
""")
    expected_contigs = ['CHROM1', 'PLASMID1']
    actual_contigs = get_vcf_contigs(fake_vcf)
    self.assertEqual(actual_contigs, expected_contigs)

if __name__ == '__main__':
  unittest.main()
