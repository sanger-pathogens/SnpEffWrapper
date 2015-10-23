#!/usr/bin/env python3

import argparse
import logging
import os
import re
import shutil
import unittest
import yaml

from jinja2 import Environment, PackageLoader

class MissingSNPEffError(ValueError):
  pass

class NoCommonContigsError(ValueError):
  pass

class MissingCodonTableError(ValueError):
  pass

class UnknownCodingTableError(ValueError):
  pass


def _java_version_ok(java):
  try:
    output = subprocess.check_output([java, '-version'], stderr=subprocess.STDOUT)
    first_line = output.decode("utf-8").splitlines()[0]
    match = re.match('^java version "1\.7\.[^"]+"$', first_line)
    return match is not None
  except CalledProcessError: # Probably this java doesn't exist
    return False
  except IndexError: # Didn't return any output
    return False

def _choose_java():
  path_java = shutil.which('java')
  if not path_java is None and _java_version_ok(path_java):
    return path_java
  sanger_pathogens_java='/software/pathogen/external/apps/usr/local/jdk1.7.0_21/bin/java'
  if not path_java is None and _java_version_ok(sanger_pathogens_java):
    return sanger_pathogens_java
  raise WrongJavaError("Could not find a suitable version of Java (1.7)")

def parse_arguments():
  # set default coding table e.g. 'default: Bacterial_and_Plant_Plastid'
  parser = argparse.ArgumentParser()
  parser.add_argument('--snpeff-exec', type=argparse.FileType('r'))
  parser.add_argument('--java-exec', type=argparse.FileType('r'))
  parser.add_argument('--coding-table', type=str,
                      default='default: Bacterial_and_Plant_Plastid')
  args = parser.parse_args()

  if args.snpeff_exec is None:
    args.snpeff_exec = shutil.which('snpeff')
  if args.snpeff_exec is None:
    raise MissingSNPEffError("Could not find snpeff in PATH")

  if args.java_exec is None:
    args.java_exec = _choose_java()
  else:
    if not _java_version_ok(args.java_exec):
      raise WrongJavaError("Needs Java 1.7, %s isn't or couldn't be found" % args.java_exec)

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

def create_temp_database(gff_file):
  database_dir = tempfile.TemporaryDirectory(prefix='snpeff_data_dir_', dir=os.getcwd())
  data_dir = os.path.join([database_dir, 'data'])
  os.mkdirs(data_dir, mode=0o644)
  shutil.copy(gff.name, os.path.join([data_dir,'genes.gff']))
  return database_dir

def get_genome_name(gff_file):
  return re.sub('\.gff(\.gz)?$', '', gff_file.name)

def create_config_file(temp_database_dir, genome_name, vcf_contigs,
                       coding_table):
  env = Environment(loader=PackageLoader('annotateVCF', 'data'))
  template = env.get_template('config.template')
  output_filename = os.path.join(temp_database_dir, 'config')
  config_content = template.render(
    temp_database_dir=temp_database_dir,
    genome_name=genome_name,
    vcf_contigs=vcf_contigs,
    coding_table=coding_table
  )
  with open(output_filename, 'w') as output_file:
    print(config_content, file=output_file, flush=True)
  return output_filename

def _snpeff_build_database(java_exec, snpeff_exec, config_filename):
  command = [java_exec, "-Xmx4g", "-jar",
             snpeff_exec, "build",
             "-gff3", "-verbose",
             "data",
             "-c", config_filename]
  logging.info("Building snpeff database")
  logging.debug("using the following command: %s" % " ".join(command))
  try:
    subprocess.check_call(command)
  except CalledProcessError:
    raise BuildDatabaseError("Problem building the database from your GFF")

def _snpeff_annotate(java_exec, snpeff_exec, vcf_filename, config_filename,
                     output_file):
  command = [java_exec, "-Xmx4g", "-jar",
             snpeff_exec, "ann",
             "-nodownload", "-verbose",
             "-c", config_filename,
             "data",
             vcf_filename]
  logging.info("Annotating %s" % vcf_filename)
  logging.debug("using the following command: %s" % " ".join(command))
  logging.debug("writing output to %s" % output_file.name)
  try:
    subprocess.check_call(command, stdout=output_file)
  except CalledProcessError:
    raise AnnotationError("Problem annotating %s" % vcf_filename)

def annotate_vcf(temp_database_dir, java_exec, snpeff_exec, vcf_file, config_filename):
  temp_output_file = tempfile.NamedTemporaryFile(mode='w', delete=False,
                                                 dir=temp_database_dir,
                                                 prefix='snpeff_annotation_',
                                                 suffix='.vcf')
  vcf_filename = vcf_file.name
  _snpeff_build_database(java_exec, snpeff_exec, config_filename)
  _snpeff_annotate(java_exec, snpeff_exec, vcf_filename, config_filename, temp_output_file)
  return temp_output_file

def annotate_vcf(args):
  coding_table = parse_coding_table(args.coding_table)
  gff_contigs = get_gff_contigs(args.gff_file)
  vcf_contigs = get_vcf_contigs(args.vcf_file)
  check_contigs(vcf_contigs, gff_contigs, coding_table)
  temp_database_dir = create_temp_database(args.gff_file)
  genome_name = get_genome_name(args.gff_file)
  config_filename = create_config_file(temp_database_dir, genome_name,
                                   vcf_contigs, coding_table)
  annotated_vcf = annotate_vcf(temp_database_dir, args.java_exec, args.snpeff_exec,
                                    args.vcf, config_filename)
  check_annotations(annotated_vcf)
  move_annotated_vcf(annotated_vcf_path, args.output_vcf.name)
  delete_temp_database(temp_database_dir)

if __name__ == '__main__':
  args = parse_arguments()
  annotate_vcf(args)
