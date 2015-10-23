#!/usr/bin/env python3

import tempfile
import unittest
import pkg_resources

from io import StringIO
from unittest.mock import patch, MagicMock

from annotateVCF.annotateVCF import *

class TestAnnotateVCF(unittest.TestCase):
  def setUp(self):
    class FakeArgs(object):
      pass
    self.fake_args = FakeArgs()

  @patch('annotateVCF.annotateVCF._java_version_ok')
  @patch('annotateVCF.annotateVCF.argparse.ArgumentParser')
  @patch('annotateVCF.annotateVCF.shutil')
  def test_snpeff_not_in_path(self, shutil_mock, argument_parser_mock, java_ok):
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

  @patch('annotateVCF.annotateVCF.logging.warn')
  def test_check_contigs(self, warn_mock):
    vcf_contigs = ['CHROM1']
    gff_contigs = ['CHROM1']
    coding_table = {'default': 'Bacterial_and_Plant_Plastid'}
    check_contigs(vcf_contigs, gff_contigs, coding_table)
    warn_mock.assert_not_called()

    vcf_contigs = ['CHROM1']
    gff_contigs = ['CHROM1']
    coding_table = {'foo': 'Bacterial_and_Plant_Plastid'}
    self.assertRaises(MissingCodonTableError, check_contigs, vcf_contigs, gff_contigs, coding_table)
    warn_mock.assert_called_once_with('Cannot annotate VCF, no coding table set for \'CHROM1\'')
    warn_mock.reset_mock()

    vcf_contigs = ['CHROM1', 'PLASMID1']
    gff_contigs = ['CHROM1', 'PLASMID1']
    coding_table = {'foo': 'Bacterial_and_Plant_Plastid'}
    self.assertRaises(MissingCodonTableError, check_contigs, vcf_contigs, gff_contigs, coding_table)
    warn_mock.assert_any_call('Cannot annotate VCF, no coding table set for \'CHROM1\'')
    warn_mock.assert_any_call('Cannot annotate VCF, no coding table set for \'PLASMID1\'')
    warn_mock.reset_mock()

    vcf_contigs = ['CHROM1', 'PLASMID1']
    gff_contigs = ['CHROM1', 'PLASMID1']
    coding_table = {'CHROM1': 'Bacterial_and_Plant_Plastid'}
    self.assertRaises(MissingCodonTableError, check_contigs, vcf_contigs, gff_contigs, coding_table)
    warn_mock.assert_any_call('Cannot annotate VCF, no coding table set for \'PLASMID1\'')
    warn_mock.reset_mock()

    vcf_contigs = ['CHROM1', 'PLASMID1']
    gff_contigs = ['ANOTHER_CHROM1', 'ANOTHER_PLASMID1']
    coding_table = {'ANOTHER_CHROM1': 'Bacterial_and_Plant_Plastid'}
    self.assertRaises(MissingCodonTableError, check_contigs, vcf_contigs, gff_contigs, coding_table)
    warn_mock.assert_any_call('Cannot annotate VCF, no coding table set for \'CHROM1\'')
    warn_mock.assert_any_call('Cannot annotate VCF, no coding table set for \'PLASMID1\'')
    warn_mock.reset_mock()

    vcf_contigs = ['CHROM1']
    gff_contigs = ['CHROM1', 'PLASMID1']
    coding_table = {'default': 'Bacterial_and_Plant_Plastid'}
    check_contigs(vcf_contigs, gff_contigs, coding_table)
    warn_mock.assert_not_called()

    vcf_contigs = ['CHROM1', 'PLASMID1']
    gff_contigs = ['CHROM1']
    coding_table = {'default': 'Bacterial_and_Plant_Plastid'}
    check_contigs(vcf_contigs, gff_contigs, coding_table)
    warn_mock.assert_any_call('Could not annotate contig \'PLASMID1\', no annotation data')
    warn_mock.reset_mock()

    vcf_contigs = ['CHROM1', 'PLASMID1', 'PLASMID2']
    gff_contigs = ['CHROM1']
    coding_table = {'default': 'Bacterial_and_Plant_Plastid'}
    check_contigs(vcf_contigs, gff_contigs, coding_table)
    warn_mock.assert_any_call('Could not annotate contig \'PLASMID1\', no annotation data')
    warn_mock.assert_any_call('Could not annotate contig \'PLASMID2\', no annotation data')
    warn_mock.reset_mock()

    vcf_contigs = ['PLASMID1', 'PLASMID2']
    gff_contigs = ['CHROM1']
    coding_table = {'default': 'Bacterial_and_Plant_Plastid'}
    self.assertRaises(NoCommonContigsError, check_contigs, vcf_contigs, gff_contigs, coding_table)
    warn_mock.assert_any_call('Could not annotate contig \'PLASMID1\', no annotation data')
    warn_mock.assert_any_call('Could not annotate contig \'PLASMID2\', no annotation data')
    warn_mock.reset_mock()

    vcf_contigs = ['CHROM1']
    gff_contigs = ['CHROM1']
    coding_table = {'default': 'UNKNOWN_CODING_TABLE'}
    self.assertRaises(UnknownCodingTableError, check_contigs, vcf_contigs, gff_contigs, coding_table)
    warn_mock.assert_any_call('Could not find coding table \'UNKNOWN_CODING_TABLE\'')
    warn_mock.reset_mock()

  def test_get_genome_name(self):
    class FakeFile(object):
      def __init__(self, name):
        self.name = name

    gff_file = FakeFile('foo')
    self.assertEqual(get_genome_name(gff_file), 'foo')

    gff_file = FakeFile('foo.gff')
    self.assertEqual(get_genome_name(gff_file), 'foo')

    gff_file = FakeFile('foo.gff.gz')
    self.assertEqual(get_genome_name(gff_file), 'foo')

    gff_file = FakeFile('foogffxgff')
    self.assertEqual(get_genome_name(gff_file), 'foogffxgff')

    gff_file = FakeFile('foogff.gff')
    self.assertEqual(get_genome_name(gff_file), 'foogff')

    gff_file = FakeFile('foo.gz')
    self.assertEqual(get_genome_name(gff_file), 'foo.gz')

  @patch('annotateVCF.annotateVCF.open', create=True)
  def test_create_config_file(self, open_mock):
    try:
      fake_output_file = tempfile.NamedTemporaryFile(mode='w',
                                                     prefix='annotateVCF_tmp_',
                                                     dir=os.getcwd(),
                                                     delete=False)
      open_mock.return_value = fake_output_file
      fake_output_filename = fake_output_file.name
  
      tests_dir = pkg_resources.resource_filename('annotateVCF', 'tests')
      expected_config_filename = os.path.join(tests_dir, 'data', 'config')
      temp_database_dir = '/tmp/fake_dir'
      genome_name = 'fake_genome'
      vcf_contigs = ['CHROM1', 'PLASMID1']
      coding_table = {
        'default': 'Standard',
        'CHROM1': 'Bacterial_and_Plant_Plastid'
      }
      config_filename = create_config_file(temp_database_dir, genome_name,
                                           vcf_contigs, coding_table)
      open_mock.called_once_with(config_filename, 'w')
      with open(fake_output_filename, 'r') as fake_output_file:
        actual_config = fake_output_file.read()
      with open(expected_config_filename, 'r') as expected_config_file:
        expected_config = expected_config_file.read()
      self.assertMultiLineEqual(actual_config, expected_config)
    except:
      os.remove(fake_output_filename)
      raise
    finally:
      os.remove(fake_output_filename)

if __name__ == '__main__':
  unittest.main()
