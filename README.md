# SnpEffWrapper
Takes a VCF and infers annotations and variant effects from a GFF using [SnpEff](http://snpeff.sourceforge.net/).

[![Build Status](https://travis-ci.org/sanger-pathogens/SnpEffWrapper.svg?branch=master)](https://travis-ci.org/sanger-pathogens/SnpEffWrapper)  
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/sanger-pathogens/SnpEffWrapper/blob/master/LICENSE)  
[![codecov](https://codecov.io/gh/sanger-pathogens/SnpEffWrapper/branch/master/graph/badge.svg)](https://codecov.io/gh/sanger-pathogens/SnpEffWrapper)

## Content
  * [Introduction](#introduction)
  * [Installation](#installation)
  * [Usage](#usage)
    * [Example usage](#example-usage)
    * [Alternative coding tables](#alternative-coding-tables)
    * [Input](#input)
  * [License](#license)
  * [Feedback/Issues](#feedbackissues)
  * [Citation](#citation)

## Introduction
SnpEff is a tool that annotates and predicts the effects of variants on genes. SnpEffWrapper takes a VCF and, using SnpEff, infers annotations and variation effects from a GFF. If you use SnpEffWrapper, please consider [citing SnpEff](http://snpeff.sourceforge.net/SnpEff.html#citing). This software is not endorsed in any respect by the original authors.

## Installation
SnpEffWrapper has the following dependencies:

 * SnpEff (>= 4.1)
 * Java (>= 1.7)
 * Jinja2
 * PyVCF
 * PyYAML

Details for the installation are provided below. If you encounter an issue when installing SnpEffWrapper please contact your local system administrator. If you encounter a bug please log it [here](https://github.com/sanger-pathogens/SnpEffWrapper/issues) or email us at path-help@sanger.ac.uk

Install [snpEff](http://snpeff.sourceforge.net/) and Java 1.7 then

```
pip install git+https://github.com/sanger-pathogens/SnpEffWrapper.git
```

## Usage
```
$ snpEffBuildAndRun --help
usage: snpEffBuildAndRun [-h] [--snpeff-exec SNPEFF_EXEC]
                         [--java-exec JAVA_EXEC] [--coding-table CODING_TABLE]
                         [-o OUTPUT_VCF] [--debug] [--keep]
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

* snpEffBuildAndRun will look for SnpEFF.jar in the following locations:
  * the file specified by `--snpeff-exec`
  * `snpEff.jar` in your local directory
  * `snpEff.jar` in your `PATH`
* SnpEff needs Java 1.7 to run; snpEffBuildAndRun will look in the following locations:
  * the file specified by `--java-exec`
  * `java` in your `PATH`

### Example usage
```
$ snpEffBuildAndRun snpEffWrapper/tests/data/minimal.gff snpEffWrapper/tests/data/minimal.vcf -o minimal.annotated.vcf
```
### Alternative coding tables

You can provide a coding table for each VCF contig otherwise it'll default to SnpEff's 'Bacterial_and_Plant_Plastid'. You can do this by providing a mapping for each contig in your VCF to the relevant table in [snpEffWrapper/data/config.template](snpEffWrapper/data/config.template) in YAML format.

For example:
```
snpEffBuildAndRun minimal.gff minimal.vcf \
  --coding-table 'default: Standard'
  
snpEffBuildAndRun minimal.gff minimal.vcf \
  --coding-table '{CHROM1: Standard, MITO1: Mitochondrial}'
  
snpEffBuildAndRun minimal.gff minimal.vcf \
  --coding-table '{default: Standard, MITO1: Mitochondrial}'
```

NB you don't need curly brackets if you're only mapping one contig (or setting a default); you do need them if you're setting different coding tables.

### Input

* The GFF must contain the reference sequence in Fasta format
* The VCF must be aligned against the reference in the GFF
* At least one of the contigs in the VCF must have annotation data in the GFF (you'll get warnings for each VCF config not in the GFF)
* You cannot provide unknown coding tables (i.e. that can't be found in [config.template](snpEffWrapper/data/config.template))

## License
SnpEffWrapper is free software, licensed under [GPLv3](https://github.com/sanger-pathogens/snpeffwrapper/blob/master/LICENSE).

## Feedback/Issues
Please report any issues to the [issues page](https://github.com/sanger-pathogens/snpeffwrapper/issues) or email path-help@sanger.ac.uk.

## Citation
If you use this, please consider [citing SnpEff](http://snpeff.sourceforge.net/SnpEff.html#citing). This software is not endorsed in any respect by the original authors.