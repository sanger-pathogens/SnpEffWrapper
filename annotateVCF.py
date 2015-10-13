#!/usr/bin/env python3

import unittest

class DatabaseCache(object):
  def __init__(self, database_dir):
    self.database_dir = database_dir
    self.create_database_dir_if_missing()
    self.create_config_yaml_if_missing()
    self.update_snpeff_config()

  def get_database(self, gff_file):
    config = self.load_config()
    database_checksum = self.get_md5(gff_file)
    try:
      database_name = config['databases'][database_checksum]
      database = AnnotationDatabase.load(database_name, self.database_dir)
    except ValueError:
      database = self._add_database(gff_file)
      config['databases'][database_checksum] = database.name
      self.save_config(config)
      self.update_snpeff_config()
    return database

  def _add_database(self, gff_file):
    database_name = self._get_database_name(gff_file)
    database = AnnotationDatabase.create(database_name, self.database_dir,
                                         gff_file)
    return database

class AnnotationDatabase(object):
  pass

if __name__ == '__main__':
  args = parse_arguments()
  database_cache = DatabaseCache(args.database_dir)
  database = database_cache.get_database(args.gff)
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
