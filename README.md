# AnnotateVCF

[![Build Status](https://travis-ci.org/sanger-pathogens/AnnotateVCF.svg?branch=master)](https://travis-ci.org/sanger-pathogens/AnnotateVCF)

Takes a VCF and applies annotations from a GFF using [SnpEff](http://snpeff.sourceforge.net/)

If you use this, please consider [citing SnpEff](http://snpeff.sourceforge.net/SnpEff.html#citing); it
made making this tool a lot easier.

## Usage

```
annotateVCF annotateVCF/tests/data/minimal.gff annotateVCF/tests/data/minimal.vcf -o minimal.annotated.vcf
```

```
$ annotateVCF --help
usage: annotateVCF [-h] [--snpeff-exec SNPEFF_EXEC] [--java-exec JAVA_EXEC]
                   [--coding-table CODING_TABLE] [-o OUTPUT_VCF] [--debug]
                   [--keep]
                   gff_file vcf_file

Takes a VCF and applies annotations from a GFF using SnpEff

positional arguments:
  gff_file              GFF with annotations including a reference genome
                        sequence
  vcf_file              VCF input to annotate (NB must be aligned to the
                        reference in your GFF

optional arguments:
  -h, --help            show this help message and exit
  --snpeff-exec SNPEFF_EXEC
                        Path to your prefered SnpEff executable (default:
                        snpEff.jar)
  --java-exec JAVA_EXEC
                        Path to Java 1.7 (default: java)
  --coding-table CODING_TABLE
                        A mapping of contig name to coding table formatted in
                        YAML
  -o OUTPUT_VCF, --output_vcf OUTPUT_VCF
                        Output for the annotated VCF (default: stdout)
  --debug               Show lots of SnpEff and other debug output
  --keep                Keep temporary files and databases (useful for
                        debugging)
```

* annotateVCF will look for SnpEFF.jar in the following locations:
  * the file specified by `--snpeff-exec`
  * `snpEff.jar` in your local directory
  * `snpEff.jar` in your `PATH`
* SnpEff needs Java 1.7 to run; annotateVCF will look in the following locations:
  * the file specified by `--java-exec`
  * `java` in your `PATH`

## Alternative coding tables

You can provide a coding table for each VCF contig else it'll default to
SnpEff's 'Bacterial_and_Plant_Plastid'.  You should provide a mapping from
the names of contigs in your VCF to the relevant table in [annotateVCF/data/config.template](annotateVCF/data/config.template)

## Built in sanity checks

* The GFF must contain the reference sequence in Fasta format
* The VCF must be aligned against the reference in the GFF
* At least one of the contigs in the VCF must have annotation data in the GFF
  (you'll get warnings for each VCF config not in the GFF)
* You cannot provide unknown coding tables (i.e. that can't be found in
  [config.template](annotateVCF/data/config.template))

## Installation

```
pip install git+https://github.com/sanger-pathogens/AnnotateVCF.git
```
