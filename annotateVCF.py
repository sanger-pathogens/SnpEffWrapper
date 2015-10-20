#!/usr/bin/env python3

import shutil
import unittest

class MissingSNPEffError(ValueError):
  pass

def parse_arguments():
  # set default coding table e.g. 'default: Bacterial_and_Plant_Plastid'
  if shutil.which('snpeff') is None:
    raise MissingSNPEffError("Could not find snpeff in PATH")

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
