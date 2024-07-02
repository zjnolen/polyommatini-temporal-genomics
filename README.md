# Polyommatini museum manuscript

This project looks at how genetic diversity has changed over time in three
species of Polyommatini butterflies, as well as a small look at a few others.
Inside this project folder is all the main working files for the data analysis
as well as the manuscript.

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
cd GenErode
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

These produce the main results files for the manuscript.

## Manuscript and Figures

The manuscript and figures were made using Quarto, primarily with R for the
figures. [`environment.yaml`](environment.yaml) is the conda environment needed
to compile the document. Notebooks for figures are stored in
[notebooks](notebooks), with the main manuscript as [`index.html`](index.html).
