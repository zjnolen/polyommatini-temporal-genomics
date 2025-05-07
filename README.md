# Species-specific loss of genetic diversity and exposure of deleterious mutations following agricultural intensification

Zachary J. Nolen, Patrycja Jamelska, Ana Sofia Torres Lara, Niklas Wahlberg,
Anna Runemark

[![Repo DOI](https://img.shields.io/badge/Repository_DOI-10.5281/zenodo.13902816-blue)](https://doi.org/10.5281/zenodo.13902816)

This repository contains the code and resources used for our study examining the
changes in genetic diversity, differentiation, inbreeding, and deleterious
mutation burdden in three species of Polyommatini butterflies in southern
Sweden. To do this, we compared genomic data from museum specimens of these
species to modern specimens, capturing a decline in genetic diversity over the
past century in a landscape characterized by agricultural intensification.

![Figure 1 from the manuscript, depicting genetic diversity decline and
inbreeding and genetic differentiation increases](figures/fig-diversity.png)

## Data Analyses

For the data analysis, the analyses are performed primarily using Snakemake,
with the main workflow in the [analyses](analyses) folder. The main workflows
uses [PopGLen](https://github.com/zjnolen/PopGLen) as a base. As the majority of
analyses performed are not natively in this workflow, it has been heavily
extended with additional Snakefiles to facilitate constructing and mapping to
variation graphs and performing genotype call based analyses. These can be found
in the [analyses/workflow/rules](analyses/workflow/rules) folder and is where
most of the code is available to look at.

By putting the analyses into Snakemake workflows,all analyses are replicable
using ~4 Snakemake runs per species - one each for constructing the variation
graphs, estimating GERP scores with GenErode, running the GONE analyses, and
running the historical/modern comparison analyses. Replication of the analyses
described in the associated manuscript must be done per species, and requires
the following steps.

### 1. Acquire data from NCBI and set up Snakemake

There are several genetic resources needed for replicating this study, described
in the supplementary material and methods text:

- The short read data of the individuals studied
- Reference genomes for the study species as well as of the outgroup species
- Functional annotations for the three focal species from Ensembl

These must be placed in the folders indicated in the config files, or config
models should be updated with their locations.

### 2. Construct variation graphs

Variation graphs were used for mapping to help reduce reference bias in
historical samples. These graphs allow the aligner to be aware of not only the
sequence of the reference genome, but also of known alternative variants. These
graphs were constructed per species using variants identified in the modern
samples. The longer read length of modern samples allows for more reliable
discovery of alternate variants, and given the short time frame we are working
with, many of these variants are likely to be present in the historical samples
and can be used to help improve their mapping.

The variation graphs are constructed with a single Snakemake command per
species:

```bash
# Run the workflow inside the analyses directory
cd analyses

snakemake --configfile config/config_<species>_make-vg.yaml \
  --snakefile workflows/Snakefile_vg-construct <additional Snakemake opts>
```

This will map the modern samples with BWA using PopGLen, call high quality
variants, then construct and index the variation graph using vg.

### 3. Estimate GERP scores

GERP scores were used to estimate the conservation of sites in the frame of the
reference genome. We used GenErode to calculate these scores, as well as infer
ancestral states to polarize the other estimates of deleterious mutation burden.
To run this, first clone the GenErode repository, checkout v0.6.2, modify the
fa2fq.py script as described in the manuscript, and then run the workflow:

```bash
git clone https://github.com/NBISweden/GenErode
# Make a folder for the species
mv GenErode generode-<species>-100
# Enter the folder to run generode
cd generode-<species>-100
# Check out the version used in the manuscript
git checkout v0.6.2
# Modify the fa2fq.py script to use 100bp fragments instead of 35
sed -i 's/35/100/g' workflow/scripts/fa2fq.py
# Run the workflow for the three focal reference genomes
snakemake --configfile ../config/generode_<reference_id>.yaml \
  <additional Snakemake opts>
```

### 4. Run GONE analyses and simulations

The GONE analyses can be run anytime after the variation graph has bee made as
they will just use the VCF used to make the graph. This is done with the the
following command:

```bash
snakemake --configfile ../config/config_<species>_vg_notrans.yaml \
  --snakefile workflows/Snakefile_gone <additional Snakemake opts>
```

The simulations with Slendr can be run after the GONE outputs are produced,
with the scripts in the [analyses/slendr](analyses/slendr) folder. See the
comments in the scripts for how they're run. I just ran these manually with
command line options and parallel per species.

### 5. Run the main analyses

Then, to map samples to the variation graph and run the main population genomic
workflow, you can run the main workflow as follows:

```bash
cd angsd
snakemake --configfile ../config/config_<species>_vg_notrans.yaml \
  <additional Snakemake opts>
```

### A few extra configs

There are a few extra configs in the config folder used for supplemental info
in the mansucript, such as running some analyses with bwa as the aligner to
compare with vg or running bwa alignments under different settings. Info about
these is in the config folder and they mostly just are run by being passed to
snakemake with the `--configfile` option.

## Manuscript and Figures

The manuscript and figures were made using Quarto, primarily with R for the
figures. [`environment.yaml`](environment.yaml) is the conda environment needed
to compile the document. If Quarto is not already installed on the machine, it
can be addded to the `environment.yaml` file. Notebooks for figures are stored
in [notebooks](notebooks).
