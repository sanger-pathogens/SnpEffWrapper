#!/usr/bin/env python3

import argparse
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
import vcf
import yaml

from collections import Counter
from jinja2 import Environment, PackageLoader
from subprocess import CalledProcessError

class MissingSNPEffError(ValueError):
  pass

class WrongJavaError(ValueError):
  pass

class NoCommonContigsError(ValueError):
  pass

class MissingCodonTableError(ValueError):
  pass

class UnknownCodingTableError(ValueError):
  pass

class BuildDatabaseError(ValueError):
  pass

class AnnotationError(ValueError):
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
    logging.debug("Using %s" % path_java)
    return path_java
  sanger_pathogens_java='/software/pathogen/external/apps/usr/local/jdk1.7.0_21/bin/java'
  if not path_java is None and _java_version_ok(sanger_pathogens_java):
    logging.debug("Using %s" % sanger_pathogens_java)
    return sanger_pathogens_java
  raise WrongJavaError("Could not find a suitable version of Java (1.7)")

def parse_arguments():
  logging.debug("Parsing user arguments")
  parser = argparse.ArgumentParser()
  parser.add_argument('--snpeff-exec', type=argparse.FileType('r'),
                     help='Path for SnpEff executable')
  parser.add_argument('--java-exec', type=argparse.FileType('r'),
                     help='Path to Java version 1.7')
  parser.add_argument('--coding-table', type=str,
                      default='default: Bacterial_and_Plant_Plastid')
  parser.add_argument('gff_file', type=argparse.FileType('r'),
                      help="GFF including reference genome sequence")
  parser.add_argument('vcf_file', type=argparse.FileType('r'),
                      help="VCF input (NB must be aligned to sequence in reference GFF")
  parser.add_argument('-o', '--output_vcf', type=argparse.FileType('w'),
                      default=sys.stdout,
                      help="Output for annotated VCF (default: stdout)")
  parser.add_argument('--verbose', action='store_false')
  args = parser.parse_args()

  if args.snpeff_exec is None:
    args.snpeff_exec = shutil.which('snpeff')
  else:
    args.snpeff_exec = args.snpeff_exec.name
  if args.snpeff_exec is None:
    raise MissingSNPEffError("Could not find snpeff in PATH")

  if args.java_exec is None:
    args.java_exec = _choose_java()
  else:
    args.java_exec = args.java_exec.name
    if not _java_version_ok(args.java_exec):
      raise WrongJavaError("Needs Java 1.7, %s isn't or couldn't be found" % args.java_exec)

  return args

def parse_coding_table(coding_table_str):
  logging.debug('Parsing the coding table')
  return yaml.load(coding_table_str)

def get_gff_contigs(gff_file):
  """Hacky gff parser to get contigs

  Just looks for the contigs, assumes they're the first column
  of a tab delimited file where the line doesn't start with '#'"""
  logging.debug('Getting the contigs from the GFF')
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
  logging.info("Checking that the VCF and GFF contigs are consistent")

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
  database_dir = tempfile.mkdtemp(prefix='snpeff_data_dir_', dir=os.getcwd())
  logging.debug("Creating directory %s for temporary database" % database_dir)
  data_dir = os.path.join(database_dir, 'data')
  logging.debug("data_dir: %s" % data_dir)
  os.makedirs(data_dir, mode=0o755)
  shutil.copy(gff_file.name, os.path.join(data_dir,'genes.gff'))
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
  logging.debug("Writing config to %s" % output_filename)
  with open(output_filename, 'w') as output_file:
    print(config_content, file=output_file, flush=True)
  return output_filename

def _snpeff_build_database(java_exec, snpeff_exec, config_filename, stdout,
                           stderr):
  command = [java_exec, "-Xmx4g", "-jar",
             snpeff_exec, "build",
             "-gff3", "-verbose",
             "data",
             "-c", config_filename]
  logging.info("Building snpeff database")
  logging.debug("using the following command: %s" % " ".join([str(c) for c in command]))
  try:
    subprocess.check_call(command, stdout=stdout, stderr=stderr)
  except CalledProcessError:
    raise BuildDatabaseError("Problem building the database from your GFF")

def _snpeff_annotate(java_exec, snpeff_exec, vcf_filename, config_filename,
                     output_file, stderr, annotation_stats_file):
  command = [java_exec, "-Xmx4g", "-jar",
             snpeff_exec, "ann",
             "-nodownload", "-verbose",
             "-stats", annotation_stats_file,
             "-c", config_filename,
             "data",
             vcf_filename]
  logging.info("Annotating %s" % vcf_filename)
  logging.debug("using the following command: %s" % " ".join(command))
  logging.debug("writing output to %s" % output_file.name)
  try:
    subprocess.check_call(command, stdout=output_file, stderr=stderr)
  except CalledProcessError:
    raise AnnotationError("Problem annotating %s" % vcf_filename)

def _get_snpeff_output_files(verbose):
  temp_output_file = tempfile.NamedTemporaryFile(mode='w', delete=False,
                                                 dir=temp_database_dir,
                                                 prefix='snpeff_output_',
                                                 suffix='.vcf')
  if verbose:
    build_stdout = tempfile.NamedTemporaryFile(mode='w',
                                                  delete=False,
                                                  dir=temp_database_dir,
                                                  prefix='snpeff_build_db_',
                                                  suffix='.o')
    build_stderr = tempfile.NamedTemporaryFile(mode='w',
                                                  delete=False,
                                                  dir=temp_database_dir,
                                                  prefix='snpeff_build_db_',
                                                  suffix='.e')
    annotate_stderr = tempfile.NamedTemporaryFile(mode='w',
                                                  delete=False,
                                                  dir=temp_database_dir,
                                                  prefix='snpeff_annotate_',
                                                  suffix='.e')
  else:
    build_stdout, build_stderr = sys.stdout, sys.stderr
    annotate_stderr = sys.stderr

  return temp_output_file, build_stdout, build_stderr, annotate_stderr

def run_snpeff(temp_database_dir, java_exec, snpeff_exec, vcf_file,
               config_filename, verbose):
  temp_output_file, build_stdout, build_stderr, annotate_stderr = _get_snpeff_output_files(verbose)
  vcf_filename = vcf_file.name
  _snpeff_build_database(java_exec, snpeff_exec, config_filename, build_stdout,
                        build_stderr)
  logging.debug("Outputting temporary VCF to %s" % vcf_filename)
  annotation_stats_file = os.path.join(temp_database_dir, 'snpEff_summary.html')
  _snpeff_annotate(java_exec, snpeff_exec, vcf_filename, config_filename,
                   temp_output_file, annotate_stderr, annotation_stats_file)
  temp_output_file.close()
  temp_output_file = open(temp_output_file.name, 'r')
  return temp_output_file

def check_annotations(annotated_vcf):
  logging.info("Checking the annotated VCF for common issues")
  error_map = {
    'WARNING_REF_DOES_NOT_MATCH_GENOME': "The reference base in your VCF didn't match the base in the GFF. Are you sure you have the right reference?",
    'WARNING_SEQUENCE_NOT_AVAILABLE': "A reference sequence was not available in your GFF. Please check that a reference sequence is available for every contig in your VCF",
    'ERROR_CHROMOSOME_NOT_FOUND': "A contig in your VCF could not be found in your GFF. Are you sure that contigs use consitent names between your input data and the reference?",
    'ERROR_OUT_OF_CHROMOSOME_RANGE': "One of your variants appears to be in a position beyond the end of the reference sequence. That's really weird, please check that you reference sequence matches your input data"
  }
  error_counter = Counter()
  annotated_vcf.seek(0)
  vcf_reader = vcf.Reader(annotated_vcf)
  for record in vcf_reader:
    annotations = ','.join(record.INFO['ANN'])
    counter_update = {error: 1 for error in error_map
                      if error in annotations}
    error_counter.update(counter_update)
  for error, count in error_counter.items():
    logging.warn("%s instances of '%s': %s" % (count, error,
                                               error_map[error]))
  if len(error_counter) > 0:
    raise AnnotationError("There were problems during the annotation, please review the warnings for details")

def move_annotated_vcf(annotated_vcf, output_vcf):
  if output_vcf is sys.stdout:
    annotated_vcf.seek(0)
    logging.debug("Writing output to stdout")
    print(annotated_vcf.read(), file=output_vcf)
  else:
    annotated_vcf.close()
    output_vcf.close()
    logging.debug("Moving annotated VCF from %s to %s" % (annotated_vcf.name,
                                                          output_vcf.name))
    shutil.move(annotated_vcf.name, output_vcf.name)

def delete_temp_database(temp_database_dir):
  logging.info("Deleting temporary files from %s" % temp_database_dir)
  shutil.rmtree(temp_database_dir)

def annotate_vcf(args):
  coding_table = parse_coding_table(args.coding_table)
  gff_contigs = get_gff_contigs(args.gff_file)
  vcf_contigs = get_vcf_contigs(args.vcf_file)
  check_contigs(vcf_contigs, gff_contigs, coding_table)
  temp_database_dir = create_temp_database(args.gff_file)
  genome_name = get_genome_name(args.gff_file)
  config_filename = create_config_file(temp_database_dir, genome_name,
                                   vcf_contigs, coding_table)
  annotated_vcf = run_snpeff(temp_database_dir, args.java_exec, args.snpeff_exec,
                             args.vcf_file, config_filename, args.verbose)
  check_annotations(annotated_vcf)
  move_annotated_vcf(annotated_vcf, args.output_vcf)
  delete_temp_database(temp_database_dir)

if __name__ == '__main__':
  args = parse_arguments()
  if args.verbose:
    logging.basicConfig(level=logging.DEBUG)
  else:
    logging.basicCongig(level=logging.INFO)
  annotate_vcf(args)
