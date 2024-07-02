# Configuration files for Polyommatini museum manuscript

## Main analyses

The config folder contains all the configuration files needed for the main
analyses. Files that are named `config_<species>.yaml` contain the configuration
for the main analyses for a given species, and set things like the reference,
what analyses to run, and with what parameters. They each are primarily to
configure the PopGLen workflow where most analyses are run, and at the top have
additional configuration settings for extensions to this analysis made by the
extra rules in the Snakemake workflow in the [../angsd](../angsd) folder. Each
of these configs also has a corresponding `samples_<species>.yaml` file which
contains the sample metadata. All use the same `units.tsv` file, which points
to the raw data paths and metadata.

## GERP Scores

GERP scores are calculated for each of the three focal species' references using
GenErode, and the config files for each are named `generode_<ref_id>.yaml`.
These GERP scores are needed for the main analyses, so the
`config_<species>.yaml` files actually reference output files from this
workflow, meaning it must be run for the references first.
