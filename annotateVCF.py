#!/usr/bin/env python3

import unittest

if __name__ == '__main__':
  args = parse_arguments()
  try:
    database = get_database_path(args.database_dir, args.gff)
  except MissingDatabase:
    database = build_database(args.database_dir, args.gff)
  database_contigs = get_database_contigs(args.database_dir, database)
  vcf_contigs = get_vcf_contigs(args.input_vcf)
  contigs_in_common = compare_contigs(database_contigs, vcf_contigs)
  if len(contigs_in_common) == 0:
    fail_missing_contigs(database_contigs, vcf_contigs)
  annotated_vcf = annotate_vcf(args.database_dir, database, args.input_vcf)
  annotation_errors = check_annotated_vcf(annotated_vcf)
  if len(annotation_errors):
    fail_annotation_errors(annotation_errors)
  output_annotated_vcf(annotated_vcf)

class TestAnnotateVCF(unittest.TestCase):
  def test_true(self):
    self.assertTrue(False)

  def test_snpeff_not_installed(self):
    pass
