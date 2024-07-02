# Calculates relative genetic load from filtered genotypes and set of positions
# with ancestral states and gerp scores. This is basically the sum of all
# derived alleles multiplied by their GERP score normalized by the total
# number of derived alleles per sample (see von Seth et al. 2021; Nat. Comm.)

library("dplyr")
library("tidyr")

# Read in CSV containing variant info and sample genotypes
df <- read.table(snakemake@input[["csv"]], header = TRUE)
colnames(df) <- c("chr", "pos", "ref", "alt", "dp", "gts")

# Read in CSV containing ancestral states and GERP scores. These have already
# been filtered to the top 1% of GERPs and to remove ancestral states of N
gerpdf <- read.table(snakemake@input[["gerp"]], header = FALSE)
colnames(gerpdf) <- c("chr", "pos", "anc", "gerp")

# Merge genotype dataframe with GERP scores
df <- merge(df, gerpdf, by = c("chr", "pos"))

# Remove positions where ancestral isn't segregating in dataset
df <- df[df$anc == df$ref | df$anc == df$alt, ]

# Read in sample list for genotypes and split the comma separated genotype
# column into per sample columns with sample IDs as names
samples <- read.table(snakemake@input[["samples"]])
colnames(samples) <- c("sample")
df <- separate(df, col = "gts", into = samples$sample, sep = ",")

# Read in population metadata for samples and subset it to samples in the
# variant CSV
sampmeta <- read.table(snakemake@input[["samplemeta"]],
    header = TRUE,
    sep = "\t"
)
sampmeta <- sampmeta[sampmeta$sample %in% samples$sample, ]
samples <- merge(samples, sampmeta, by = "sample", sort = FALSE)

# Convert genotypes to counts of alternate alleles
df[df == "0/0"] <- 0
df[df == "0/1"] <- 1
df[df == "1/1"] <- 2
df[df == "./."] <- NA
df <- df %>% mutate_at(
    colnames(df)[!colnames(df) %in% c(
        "chr", "pos", "ref", "alt", "dp", "gts", "anc", "gerp"
    )],
    as.numeric
)

# In situations where the alternate == the ancestral, swap the genotypes (i.e.
# 0/0 becomes 1/1 and vice versa). This converts the 0, 1, 2 genotypes to
# derived allele counts.
df[df$alt == df$anc, samples$sample] <- abs(
    df[df$alt == df$anc, samples$sample] - 2
)

# estimate relative genetic load per sample
relload <- c()

for (s in samples$sample) {
    samdf <- df[, c(s, "gerp")]
    samdf <- samdf[complete.cases(samdf), ]
    # load is per sample the total sum of the derived allele count at each
    # position * the gerp score for the position, divided by sum of the total
    # derived alleles. Excludes positions per samples where sample is missing
    # data
    sum_corr_gerp <- sum(samdf[, 1] * samdf$gerp)
    sum_der <- sum(samdf[, 1])
    load <- sum_corr_gerp / sum_der
    nsnps <- nrow(samdf)
    row <- c(s, sum_corr_gerp, sum_der, load, nsnps)
    relload <- rbind(relload, row)
}

# Add in the sample metadata to get a useful table for plotting
colnames(relload) <- c(
    "sample", "sum_corrected_gerp", "sum_derived", "relload", "nsnps"
)
relload <- merge(samples, relload, by = "sample")

# write table to output file
write.table(relload,
    file = snakemake@output[["load"]], quote = FALSE,
    sep = "\t", row.names = FALSE, col.names = TRUE
)

# Now, we generate some bootstraps for this by randomly sampling half the SNPs.
# This will help account for variation across genome and different sequencing
# quality (i.e. number of usable SNPs) between the samples.

nboot <- 300 # how many bootstraps to perform
bootsize <- as.integer(nrow(df) / 3) # how many SNPs to sample per bootstrap
bootload <- c()

for (s in samples$sample) {
    for (b in c(1:nboot)) {
        samdf <- df[, c(s, "gerp")]
        samdf <- samdf[complete.cases(samdf), ]
        samdf <- samdf[sample(nrow(samdf), bootsize), ]
        sum_corr_gerp <- sum(samdf[, 1] * samdf$gerp)
        sum_der <- sum(samdf[, 1])
        load <- sum_corr_gerp / sum_der
        nsnps <- nrow(samdf)
        row <- c(s, b, sum_corr_gerp, sum_der, load, nsnps)
        bootload <- rbind(bootload, row)
    }
}

# Add in the sample metadata to get a useful table for plotting
colnames(bootload) <- c(
    "sample", "bootstrap", "sum_corrected_gerp",
    "sum_derived", "relload", "nsnps"
)
bootload <- merge(samples, bootload, by = "sample")

# write table to output file
write.table(bootload,
    file = snakemake@output[["load_boot"]], quote = FALSE,
    sep = "\t", row.names = FALSE, col.names = TRUE
)

# As an extra test, we can use only sites with no missing data.

df <- df[complete.cases(df), ]

relload <- c()

for (s in samples$sample) {
    samdf <- df[, c(s, "gerp")]
    samdf <- samdf[complete.cases(samdf), ]
    # load is per sample the total sum of the derived allele count at each
    # position * the gerp score for the position, divided by sum of the total
    # derived alleles. Excludes positions per samples where sample is missing
    # data
    sum_corr_gerp <- sum(samdf[, 1] * samdf$gerp)
    sum_der <- sum(samdf[, 1])
    load <- sum_corr_gerp / sum_der
    nsnps <- nrow(samdf)
    row <- c(s, sum_corr_gerp, sum_der, load, nsnps)
    relload <- rbind(relload, row)
}

# Add in the sample metadata to get a useful table for plotting
colnames(relload) <- c(
    "sample", "sum_corrected_gerp", "sum_derived", "relload", "nsnps"
)
relload <- merge(samples, relload, by = "sample")

# write table to output file
write.table(relload,
    file = snakemake@output[["load_nomiss"]], quote = FALSE,
    sep = "\t", row.names = FALSE, col.names = TRUE
)
