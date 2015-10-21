#!/usr/bin/env python3

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

if __name__ == '__main__':
  unittest.main()
