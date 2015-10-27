#!/usr/bin/env python3

import tempfile
import unittest
import pkg_resources
import vcf

from io import StringIO
from unittest.mock import patch, MagicMock

from annotateVCF.annotateVCF import *
from annotateVCF.annotateVCF import _java_version_ok, _choose_java

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
    parsed_args.snpeff_exec.name = 'foobar'
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

  @patch('annotateVCF.annotateVCF.shutil.which')
  @patch('annotateVCF.annotateVCF._java_version_ok')
  def test_choose_java(self, java_ok_mock, which_mock):
    which_mock.return_value = '/foo/bar/java'
    java_ok_mock.side_effect = lambda java: java in list_of_ok_javas

    list_of_ok_javas = [
      '/foo/bar/java'
    ]
    self.assertEqual(_choose_java(), '/foo/bar/java')

    list_of_ok_javas = [
      '/software/pathogen/external/apps/usr/local/jdk1.7.0_21/bin/java'
    ]
    self.assertEqual(_choose_java(),
                     '/software/pathogen/external/apps/usr/local/jdk1.7.0_21/bin/java')

    list_of_ok_javas = [
      '/foo/bar/java',
      '/software/pathogen/external/apps/usr/local/jdk1.7.0_21/bin/java'
    ]
    self.assertEqual(_choose_java(), '/foo/bar/java')

    list_of_ok_javas = []
    self.assertRaises(WrongJavaError, _choose_java)

  @patch('annotateVCF.annotateVCF.logging.warn')
  def test_check_annotations(self, warn_mock):
    fake_vcf = StringIO("""\
##fileformat=VCFv4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample_1	sample_2
CHROM1	400	.	G	A	.	.	ANN=A|foo|bar|	GT	0	1
""")
    self.assertEqual(check_annotations(fake_vcf), None)
    warn_mock.assert_not_called()

    fake_vcf = StringIO("""\
##fileformat=VCFv4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample_1	sample_2
CHROM1	400	.	G	A	.	.	ANN=A|foo|bar|ERROR_CHROMOSOME_NOT_FOUND	GT	0	1
""")
    self.assertRaises(AnnotationError, check_annotations, fake_vcf)
    expected_warning = "1 instances of 'ERROR_CHROMOSOME_NOT_FOUND': A contig in your VCF could not be found in your GFF. Are you sure that contigs use consitent names between your input data and the reference?"
    warn_mock.assert_called_once_with(expected_warning)

  @patch('annotateVCF.annotateVCF.delete_temp_database')
  @patch('annotateVCF.annotateVCF._get_snpeff_output_files')
  def test_happy_case(self, output_mock, delete_database_mock):
    delete_database_mock.side_effect = shutil.rmtree
    temp_annotated_vcf = tempfile.NamedTemporaryFile(mode='w', delete=False,
                                                     dir=os.getcwd(),
                                                     prefix='temp_annotated_',
                                                     suffix='.vcf')
    other_temp_files = [tempfile.NamedTemporaryFile(mode='w', delete=False,
                                                    dir=os.getcwd(),
                                                    prefix='temp_output_')
                        for i in range(3)]
    output_mock.return_value = [temp_annotated_vcf] + other_temp_files

    fake_args = MagicMock()
    snpeff_exec = os.environ.get('SNPEFF_EXEC', 'snpEff.jar')
    self.assertTrue(os.path.isfile(snpeff_exec), "Couldn't find snpeff, set with SNPEFF_EXEC envionment variable")
    fake_args.snpeff_exec = snpeff_exec

    java_exec = os.environ.get('JAVA_EXEC', shutil.which('java'))
    self.assertTrue(os.path.isfile(java_exec), "Could not find Java, set with JAVA_EXEC environment variable")
    self.assertTrue(_java_version_ok(java_exec), "%s isn't Java 1.7, set with JAVA_EXEC environment variable" % java_exec)
    fake_args.java_exec = java_exec

    fake_args.coding_table = 'default: Bacterial_and_Plant_Plastid'

    tests_dir = pkg_resources.resource_filename('annotateVCF', 'tests')
    minimal_gff_filename = os.path.join(tests_dir, 'data', 'minimal.gff')
    minimal_vcf_filename = os.path.join(tests_dir, 'data', 'minimal.vcf')
    fake_args.gff_file = open(minimal_gff_filename, 'r')
    fake_args.vcf_file = open(minimal_vcf_filename, 'r')

    output_annotated_vcf = tempfile.NamedTemporaryFile(mode='w', delete=False,
                                                     dir=os.getcwd(),
                                                     prefix='output_annotated_',
                                                     suffix='.vcf')
    fake_args.output_vcf = output_annotated_vcf
    fake_args.verbose = False

    annotate_vcf(fake_args)

    output_annotated_vcf = open(output_annotated_vcf.name, 'r')
    vcf_reader = vcf.Reader(output_annotated_vcf)
    (only_record,) = vcf_reader # Will raise an exception if there isn't exactly one
    (only_annotation,) = only_record.INFO['ANN'] # likewise
    self.assertRegex(only_annotation,
                     '^A\|missense_variant\|.+\|c\.14C>A\|p\.Pro5His\|')

    self.assertEqual(delete_database_mock.call_count, 1)

    for f in other_temp_files:
      f.close()
      os.remove(f.name)
    output_annotated_vcf.close()
    os.remove(output_annotated_vcf.name)

if __name__ == '__main__':
  unittest.main()
