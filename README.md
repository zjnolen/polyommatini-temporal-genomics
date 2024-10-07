# Species-specific loss of genetic diversity and accumulation of genetic load following agricultural intensification

Zachary J. Nolen, Patrycja Jamelska, Ana Sofia Torres Lara, Niklas Wahlberg,
Anna Runemark

This repository contains the code and resources used for our study examining the
changes in genetic diversity, differentiation, inbreeding, and genetic load in
three species of Polyommatini butterflies in southern Sweden. To do this, we
compared genomic data from museum specimens of these species to modern
specimens, capturing a decline in genetic diversity over the past century in a
landscape characterized by agricultural intensification.

## Data Analyses

For the data analysis, the analyses are performed primarily using Snakemake,
with the main workflow in the [angsd](angsd) folder. The main workflow uses
[PopGLen](https://github.com/zjnolen/PopGLen) as its base, and via additional
Snakefiles, extends it to perform the remaining analyses not included in the
PopGLen workflow. One of these extended analyses uses GERP scores for each
reference genome, which we calculated using several Lepidoptera references and
the [GenErode pipeline](https://github.com/NBISweden/GenErode). As the complete
manuscript workflow in the angsd folder requires the output of GenErode, the
GERP scores should be calculated with GenErode first when reproducing the
analyses.

Configurations for both workflows are stored in the [config](config) folder.
To run the GenErode workflow, clone it in the main working directory and run it
for a reference with Snakemake, using the `--configfile` option, as shown:

```bash
git clone https://github.com/NBISweden/GenErode
mv GenErode generode-<species>
cd generode-<species>
git checkout v0.6.0
snakemake --configfile ../config/generode_<reference_id>.yaml
```

Then, to run the remaining analyses for that species, run the workflow in the
angsd folder using the same `--configfile` option, pointing to the config file
required for that workflow:

```bash
cd angsd
snakemake --configfile ../config/config_<species>.yaml
```

These produce the main results files for the manuscript. A subset of these
outputs will be included in this repository for visualization purposes, size
permitting.

## Manuscript and Figures

The manuscript and figures were made using Quarto, primarily with R for the
figures. [`environment.yaml`](environment.yaml) is the conda environment needed
to compile the document. If Quarto is not already installed on the machine, it
can be addded to the `environment.yaml` file. Notebooks for figures are stored
in [notebooks](notebooks).
