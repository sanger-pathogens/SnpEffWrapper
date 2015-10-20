#!/usr/bin/env python3

import unittest

class DatabaseCache(object):
  """Representation of all of the SNPEff Databases

  Handles things like config updates and making sure we
  don't go creating duplicate databases for the same
  GFF file"""
  def __init__(self, database_dir):
    """Object representing SNPEff databases"""
    self.database_dir = database_dir
    self._create_database_dir_if_missing()
    self._create_config_yaml_if_missing()
    self._update_snpeff_config()

  def get_database(self, gff_file):
    """Get a SNPEff database if it exists

    Return an error if no matching database can be found"""
    config = self.load_config()
    database_checksum = self.get_md5(gff_file)
    try:
      database_name = config['databases'][database_checksum]
      database = AnnotationDatabase(database_name, self.database_dir)
    except ValueError:
      raise MissingDatabaseError("Could not find database for %s in cache" %
                                 gff_file.name)
    return database

  def add_database(self, gff_file, coding_table_string=None):
    """Take a GFF File and optional coding table and create a SNPEff database

    Checks if the database already exists and creates a new
    one if it doesn't"""
    database_name = self._get_database_name(gff_file)
    database_checksum = self.get_md5(gff_file)
    if database_checksum in config['databases']:
      raise DatabaseExistsError("A database for %s already exists" %
                                database_name)
    config['databases'][database_checksum] = database_name
    # copy {gff_file} to {database_dir}/{dabase_name}/genes.gff
    # parse out the names of the contigs and save them to a file
    # parse the codon table
    coding_table = self._parse_coding_table(coding_table_string)
    # update the config
    self._save_config(config)
    self._update_snpeff_config(config)
    # build the command
    # run the command
    # check the output
    return AnnotationDatabase(database_name, self.database_dir)

  def _parse_coding_table(self, coding_table_string):
    """Parses a coding table string into a map of contig to coding table

    Must either have a coding table for each contig or provide a default
    which is used for contigs without a coding table.  Coding tables
    must be 'known'"""
    pass

  def _update_snpeff_config(self, config):
    # Set the coding table for each config (using the default if present)
    pass

class AnnotationDatabase(object):
  def __init__(self, database_name, database_dir):
    pass

  def check_contigs(self, input_vcf):
    # Check that all of the contigs in the VCF are in the GFF
    pass

  def annotate_vcf(self, input_vcf):
    pass

class AnnotatedVcf(object):
  def __init__(self, file_path):
    pass

  def move(self, dest_path):
    pass

  def check_annotations(self):
    pass

if __name__ == '__main__':
  args = parse_arguments()
  # set default coding table e.g. '{"default": "Bacterial_and_Plant_Plastid"}'
  database_cache = DatabaseCache(args.database_dir)
  try:
    database = database_cache.get_database(args.gff)
  except MissingDatabaseError:
    database = database_cache.add_database(args.gff, args.coding_table)
  contigs_in_common = database.check_contigs(args.input_vcf)
  annotated_vcf = database.annotate_vcf(args.input_vcf)
  annotation_errors = annotated_vcf.check_annotations()
  if annotation_errors:
    fail_annotation_errors(annotation_errors)
  annotated_vcf.move(args.output_vcf.name)

class TestAnnotateVCF(unittest.TestCase):
  def test_true(self):
    self.assertTrue(False)

  def test_snpeff_not_installed(self):
    pass
