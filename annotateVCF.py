#!/usr/bin/env python3

import unittest

def parse_arguments():
  # set default coding table e.g. 'default: Bacterial_and_Plant_Plastid'
  pass

def annotate_vcf(args):
  coding_table = parse_coding_table(args.coding_table)
  gff_contigs = get_gff_configs(args.gff_file)
  vcf_contigs = get_vcf_configs(args.vcf_file)
  check_contigs(vcf_configs, gff_configs, coding_table)
  temp_database_dir = create_temp_database(args.data_dir, args.gff_file)
  config_file = create_config_file(temp_database_dir, args.gff_file,
                                   vcf_contigs, coding_table)
  annotated_vcf_path = annotate_vcf(args.vcf, config_file)
  check_annotations(annotated_vcf_path)
  move_annotated_vcf(annotated_vcf_path, args.output_vcf.name)
  delete_temp_database(temp_database_dir)

if __name__ == '__main__':
  args = parse_arguments()
  annotate_vcf(args)

class TestAnnotateVCF(unittest.TestCase):
  def test_true(self):
    self.assertTrue(False)

  def test_snpeff_not_installed(self):
    pass
