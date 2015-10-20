#!/usr/bin/env python3

from unittest.mock import patch

from annotateVCF import *

class TestAnnotateVCF(unittest.TestCase):
  def setUp(self):
    class FakeArgs(object):
      pass
    self.fake_args = FakeArgs()

  @patch('annotateVCF.shutil')
  def test_snpeff_not_installed(self, shutil_mock):
    shutil_mock.which.return_value = None
    self.assertRaises(MissingSNPEffError, parse_arguments)

if __name__ == '__main__':
  unittest.main()
