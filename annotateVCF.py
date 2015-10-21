#!/usr/bin/env python3

import argparse
import logging
import shutil
import unittest
import yaml

class MissingSNPEffError(ValueError):
  pass

class NoCommonContigsError(ValueError):
  pass

class MissingCodonTableError(ValueError):
  pass

class UnknownCodingTableError(ValueError):
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

def get_vcf_contigs(vcf_file):
  """Hacky vcf parser to get contigs

  Just looks for the contigs, assumes they're the first column
  of a tab delimited file where the line doesn't start with '#'"""
  vcf_file.seek(0)
  contigs = set()
  for line in vcf_file:
    if line[0] == '#':
      continue
    contig = line.split('\t')[0].strip()
    contigs.add(contig)
  return sorted(contigs)

def check_contigs(vcf_contigs, gff_contigs, coding_table):
  """Check that contigs are consistent

  If any contig in the VCF isn't in the coding table, fail.
  If not all of the VCF contigs are in the GFF, raise warnings.
  If none of the VCF contigs are in the GFF, fail"""

  # Check the VCF contigs are consistent with the coding table
  missing_coding_tables = []
  if 'default' not in coding_table:
    missing_coding_tables = [table for table in vcf_contigs
                             if table not in coding_table]
  for table in missing_coding_tables:
    logging.warn("Cannot annotate VCF, no coding table set for '%s'" % table)

  # Check the VCF contigs are consistent with the GFF contigs
  missing_contigs = [contig for contig in vcf_contigs
                     if contig not in gff_contigs]
  for contig in missing_contigs:
    logging.warn("Could not annotate contig '%s', no annotation data" % contig)

  # Check the coding_table has known encodings
  known_encodings = [
    'Alternative_Flatworm_Mitochondrial',
    'Alternative_Yeast_Nuclear',
    'Ascidian_Mitochondrial',
    'Bacterial_and_Plant_Plastid',
    'Blepharisma_Macronuclear',
    'Chlorophycean_Mitochondrial',
    'Ciliate_Nuclear',
    'Coelenterate',
    'Dasycladacean_Nuclear',
    'Echinoderm_Mitochondrial',
    'Euplotid_Nuclear',
    'Flatworm_Mitochondrial',
    'Hexamita_Nuclear',
    'Invertebrate_Mitochondrial',
    'Mitochondrial',
    'Mold_Mitochondrial',
    'Mycoplasma',
    'Protozoan_Mitochondrial',
    'Scenedesmus_obliquus_Mitochondrial',
    'Spiroplasma',
    'Standard',
    'Thraustochytrium_Mitochondrial',
    'Trematode_Mitochondrial',
    'Vertebrate_Mitochondrial',
    'Yeast_Mitochondrial'
  ]
  unknown_encodings = [enc for enc in coding_table.values()
                       if enc not in known_encodings]
  for encoding in unknown_encodings:
    logging.warn("Could not find coding table '%s'" %
                 encoding)

  # Blow up for critical issues
  if len(missing_coding_tables) > 0:
    raise MissingCodonTableError("Could not find coding tables for all contigs, see warnings for details")
  if missing_contigs == vcf_contigs:
    raise NoCommonContigsError("Could not find anotation data for any contigs, see warnings for details")
  if len(unknown_encodings) > 0:
    raise UnknownCodingTableError("Could not find coding table, see warnings for details")

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
