#!/usr/bin/env python3

import argparse
import shutil
import unittest
import yaml

class MissingSNPEffError(ValueError):
  pass

def parse_arguments():
  # set default coding table e.g. 'default: Bacterial_and_Plant_Plastid'
  parser = argparse.ArgumentParser()
  parser.add_argument('--snpeff-exec', type=argparse.FileType('r'))
  parser.add_argument('--coding-table', type=str,
                      default='default: Bacterial_and_Plant_Plastid')
  args = parser.parse_args()

  if args.snpeff_exec is None:
    args.snpeff_exec = shutil.which('snpeff')
  if args.snpeff_exec is None:
    raise MissingSNPEffError("Could not find snpeff in PATH")

  return args

def parse_coding_table(coding_table_str):
  return yaml.load(coding_table_str)

def get_gff_contigs(gff_file):
  """Hacky gff parser to get contigs

  Just looks for the contigs, assumes they're the first column
  of a tab delimited file where the line doesn't start with '#'"""
  gff_file.seek(0)
  contigs = set()
  for line in gff_file:
    if line[0] == '#':
      continue
    contig = line.split('\t')[0].strip()
    contigs.add(contig)
  return sorted(contigs)

def annotate_vcf(args):
  coding_table = parse_coding_table(args.coding_table)
  gff_contigs = get_gff_contigs(args.gff_file)
  vcf_contigs = get_vcf_contigs(args.vcf_file)
  check_contigs(vcf_contigs, gff_contigs, coding_table)
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
