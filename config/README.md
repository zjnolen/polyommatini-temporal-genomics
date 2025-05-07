# Configuration files for Polyommatini museum manuscript

## Configuration files

There are many configurations for the different parts of the analyses. These are
the formats of the filenames that you can expect in here:

- `generode_<refname>.yaml` contains the configuration for the GenErode pipeline
  to estimate the GERP scores and ancestral states.

- `config_<species>_make-vg.yaml` contains the configuration for constructing
  the variation graph from bwa mem aligned modern samples.

- `config_<species>_vg.yaml`/`config_<species>_vg_notrans.yaml` contains the
  configuration for running the main manuscript analyses. These align the
  samples to the variation graph and then run the population genomic analyses.

- `config_<species>_bwa_notrans.yaml` is the same as the previous, but aligns
  the samples with bwa mem/aln instead. This was used for comparisons of bwa and
  vg in the supplementary material.

- `config_picarus_bwa_notrans-aln<#>.yaml` is used to map a *Po. icarus*
  historical sample to the reference with alternate bwa aln settings for
  comparisons shown in the supplement.

## Sample files

Sample files contain the sample lists for each configuration and are referenced
in the respective config.yaml files

## Units files

Units files point to the raw data for the samples. These are the raw reads in
the case of `units.tsv` and are used to get the raw reads whenever they are
mapped. `units-vgmap.tsv` contains the same info, but also an added bam file so
that the workflow will know to use the vg aligned bams instead of aligning the
reads with bwa.

## GONE files

These are templates for the GONE input param files required for the GONE
analyses for each species and follow the format expected in that tool.
