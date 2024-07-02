# Polyommatini museum manuscript

This project looks at how genetic diversity has changed over time in three
species of Polyommatini butterflies, as well as a small look at a few others.
Inside this project folder is all the main working files for the data analysis
as well as the manuscript.

## Data Analyses

For the data analysis, the analysis is performed primarily using Snakemake, with
the main workflow in the [angsd](angsd) folder. This is because the main
workflow uses [PopGLen](https://github.com/zjnolen/PopGLen), our pipeline for GL
based pop gen analyses, which we extend with additional Snakemake rules. Before
running that pipeline, we generated GERP scores for each of the focal species'
reference genomes using several Lepidoptera references and the
[GenErode pipeline](https://github.com/NBISweden/GenErode). The code for that
process is available in the [generode](generode) folder, which is a clone of
GenErode at v0.6.0. Configurations for both workflows are stored in the
[config](config) folder here. To run the workflow for each, enter either the
[angsd](angsd) or [generode](generode) folder and use snakemake with the
`--configfile` option, as shown:

```bash
# For the angsd folder:
snakemake --configfile ../config/config_<species>.yaml

# For the generode folder:
snakemake --configfile ../config/generode_<reference_id>.yaml
```

## Manuscript and Figures

The manuscript and figures were made using Quarto, primarily with R for the
figures. [`environment.yaml`](environment.yaml) is the conda environment needed
to compile the document. Notebooks for figures are stored in
[notebooks](notebooks), with the main manuscript as [`index.html`](index.html).
