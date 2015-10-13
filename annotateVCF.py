#!/usr/bin/env python3

import unittest

class DatabaseCache(object):
  def __init__(self, database_dir):
    pass

  def get_database(self, gff_file):
    pass

  def add_database(self, gff_file):
    pass

class AnnotationDatabase(object):
  pass

if __name__ == '__main__':
  args = parse_arguments()
  database_cache = DatabaseCache(args.database_dir)
  try:
    database = database_cache.get_database(args.gff)
  except MissingDatabase:
    database = database_cache.add_database(args.gff)
  contigs_in_common = database.check_contigs(args.input_vcf)
  annotated_vcf = database.annotate_vcf(args.input_vcf)
  annotation_errors = check_annotated_vcf(annotated_vcf)
  if annotation_errors:
    fail_annotation_errors(annotation_errors)
  output_annotated_vcf(annotated_vcf)

class TestAnnotateVCF(unittest.TestCase):
  def test_true(self):
    self.assertTrue(False)

  def test_snpeff_not_installed(self):
    pass
