library(dplyr)
library(tidyr)

# Read in VEP output for variants
vep <- read.table(snakemake@input[["vep"]], header = FALSE)
colnames(vep) <- c("chr", "pos", "vepalt", "consequence", "notes")
vep <- separate(vep,
  col = "notes", into = c("impact"), extra = "drop",
  sep = ";"
)

# Read in variants
vars <- read.table(snakemake@input[["vars"]], header = TRUE)
colnames(vars) <- c("chr", "pos", "ref", "alt", "dp", "gt")

# Merge VEP impacts with variants and clean up impact column
df <- merge(vars, vep, by = c("chr", "pos"))
sum(df$alt == df$vepalt) == nrow(df) # Check that alt alleles match
df$impact <- gsub("IMPACT=", "", df$impact)

# Remove modifier variants
df <- df[df$impact != "MODIFIER", ]

# Get sample list and metadata
samples <- read.table(snakemake@input[["samples"]])
colnames(samples) <- c("sample")
sampmeta <- read.table(snakemake@input[["pops"]],
  header = TRUE, fill = TRUE,
  sep = "\t", comment.char = "#"
)
sampmeta <- sampmeta[, c(1:4)]
sampmeta <- sampmeta[sampmeta$sample %in% samples$sample, ]
sampord <- samples$sample
samples <- merge(samples, sampmeta, by = "sample", sort = FALSE)
samples$sample == sampord

# Split genotypes into sample columns
df <- separate(df, col = "gt", into = samples$sample, sep = ",")
# convert genotypes to alternate allele counts
df[df == "0/0"] <- 0
df[df == "0/1"] <- 1
df[df == "1/1"] <- 2
df[df == "./."] <- NA
df <- df %>%
  mutate_at(
    colnames(df)[!colnames(df) %in% c(
      "chr", "pos", "ref", "alt", "dp", "vepalt", "consequence", "impact",
      "anc", "gerp"
    )],
    as.numeric
  )

# make data frames of historical and modern genotypes only, as well as list of
# variant impacts
histgt <- df[, which(
  names(df) %in% samples$sample[sampmeta$time == "historical"]
)]
modgt <- df[, which(names(df) %in% samples$sample[sampmeta$time == "modern"])]
impacts <- df$impact

# Calculate Rxy for the dataset, treating historical samples and modern samples
# as single units

# set values for 'bootstraps'
nboot <- 1000 # number of 'bootstraps' to perform
nsample <- 8 # number of samples to sample per group

# initialize rxy data frame
rxy <- data.frame()

for (b in c(1:nboot)) {
  # subset historical genotypes by nsample
  subhistgt <- sample(histgt, nsample)
  # calculate allele frequencies for historical samples at all sites
  hist_af <- (rowSums(subhistgt, na.rm = TRUE)) /
    (rowSums(!is.na(subhistgt)) * 2)
  # subset modern genotypes by nsample
  submodgt <- sample(modgt, nsample)
  # calculate allele frequencies for modern samples at all sites
  mod_af <- (rowSums(submodgt, na.rm = TRUE)) /
    (rowSums(!is.na(submodgt)) * 2)
  # remove sites with missing individuals from subsampling to ensure subsamples
  # are even (only relevant if data wasn't already filtered for missingness)
  sitefilt <- (complete.cases(submodgt) * complete.cases(subhistgt)) == 1
  # Filter by missingness. Also only use sites that are segregating in the
  # subsample. This makes the subsample a filtered set with equivalent
  # processing to our original data
  sitefilt <- sitefilt & (
    (mod_af > 0 & mod_af < 1) | (hist_af > 0 & hist_af < 1)
  )
  hist_af <- hist_af[sitefilt]
  mod_af <- mod_af[sitefilt]
  bootimpacts <- impacts[sitefilt]
  usable_sites <- length(bootimpacts)
  # Prep Rxy dataframe
  rxy_df <- data.frame(matrix(nrow = length(bootimpacts), ncol = 0))
  rxy_df$hist_af <- hist_af
  rxy_df$mod_af <- mod_af
  rxy_df$impact <- bootimpacts
  usable_sites <- nrow(rxy_df)
  # calculate Rxy for historical-modern frequencies (<1 is increase in freq
  # over time)
  rxy_df <- rxy_df %>%
    group_by(impact) %>%
    summarise(
      rxy = ((sum(hist_af * (1 - mod_af))) / (sum(mod_af * (1 - hist_af))))
    )
  rxy_df$boot <- b
  rxy_df$usablesites <- usable_sites
  rxy_df$nsamples <- nsample
  rxy <- rbind(rxy, rxy_df)
}

# Column names update
colnames(rxy) <- c(
  "var.impact", "rxy.hist.mod", "boot", "usable.sites", "nsamples"
)

write.table(rxy,
  file = snakemake@output[["rxy"]], quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)

# Now count the totals per individual for alternate alleles. Also count the
# total alt count per individual for normalization into 'relative' load

varcounts <- data.frame()

for (sam in samples$sample) {
  high <- sum(df[df$impact == "HIGH", sam], na.rm = TRUE)
  high_hom <- sum(df[df$impact == "HIGH", sam] == 2, na.rm = TRUE) * 2
  mod <- sum(df[df$impact == "MODERATE", sam], na.rm = TRUE)
  mod_hom <- sum(df[df$impact == "MODERATE", sam] == 2, na.rm = TRUE) * 2
  low <- sum(df[df$impact == "LOW", sam], na.rm = TRUE)
  low_hom <- sum(df[df$impact == "LOW", sam] == 2, na.rm = TRUE) * 2
  alts <- sum(df[, sam], na.rm = TRUE)
  nsites <- sum(!is.na(df[, sam]), na.rm = TRUE)
  row <- c(sam, high, high_hom, mod, mod_hom, low, low_hom, alts, nsites)
  varcounts <- rbind(varcounts, row)
}

colnames(varcounts) <- c(
  "sample", "high", "high_hom", "mod", "mod_hom", "low", "low_hom", "alts",
  "nsites"
)
numcols <- c(
  "high", "high_hom", "mod", "mod_hom", "low", "low_hom", "alts",
  "nsites"
)
varcounts[numcols] <- sapply(varcounts[numcols], as.numeric)
varcounts <- merge(varcounts, samples, by = "sample")

write.table(varcounts,
  file = snakemake@output[["varcounts"]], quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)

# Limit to only sites with data for all individuals and repeat both
df <- df[complete.cases(df),]

# make data frames of historical and modern genotypes only, as well as list of
# variant impacts
histgt <- df[, which(
  names(df) %in% samples$sample[sampmeta$time == "historical"]
)]
modgt <- df[, which(names(df) %in% samples$sample[sampmeta$time == "modern"])]
impacts <- df$impact

# Calculate Rxy for the dataset, treating historical samples and modern samples
# as single units

# set values for 'bootstraps'
nboot <- 1000 # number of 'bootstraps' to perform
nsample <- 8 # number of samples to sample per group

# initialize rxy data frame
rxy <- data.frame()

for (b in c(1:nboot)) {
  # subset historical genotypes by nsample
  subhistgt <- sample(histgt, nsample)
  # calculate allele frequencies for historical samples at all sites
  hist_af <- (rowSums(subhistgt, na.rm = TRUE)) /
    (rowSums(!is.na(subhistgt)) * 2)
  # subset modern genotypes by nsample
  submodgt <- sample(modgt, nsample)
  # calculate allele frequencies for modern samples at all sites
  mod_af <- (rowSums(submodgt, na.rm = TRUE)) /
    (rowSums(!is.na(submodgt)) * 2)
  # remove sites with missing individuals from subsampling to ensure subsamples
  # are even (only relevant if data wasn't already filtered for missingness)
  sitefilt <- (complete.cases(submodgt) * complete.cases(subhistgt)) == 1
  # Filter by missingness. Also only use sites that are segregating in the
  # subsample. This makes the subsample a filtered set with equivalent
  # processing to our original data
  sitefilt <- sitefilt & (
    (mod_af > 0 & mod_af < 1) | (hist_af > 0 & hist_af < 1)
  )
  hist_af <- hist_af[sitefilt]
  mod_af <- mod_af[sitefilt]
  bootimpacts <- impacts[sitefilt]
  usable_sites <- length(bootimpacts)
  # Prep Rxy dataframe
  rxy_df <- data.frame(matrix(nrow = length(bootimpacts), ncol = 0))
  rxy_df$hist_af <- hist_af
  rxy_df$mod_af <- mod_af
  rxy_df$impact <- bootimpacts
  usable_sites <- nrow(rxy_df)
  # calculate Rxy for historical-modern frequencies (<1 is increase in freq
  # over time)
  rxy_df <- rxy_df %>%
    group_by(impact) %>%
    summarise(
      rxy = ((sum(hist_af * (1 - mod_af))) / (sum(mod_af * (1 - hist_af))))
    )
  rxy_df$boot <- b
  rxy_df$usablesites <- usable_sites
  rxy_df$nsamples <- nsample
  rxy <- rbind(rxy, rxy_df)
}

# Column names update
colnames(rxy) <- c(
  "var.impact", "rxy.hist.mod", "boot", "usable.sites", "nsamples"
)

write.table(rxy,
  file = snakemake@output[["rxy_nomiss"]], quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)

# Now count the totals per individual for alternate alleles. Also count the
# total alt count per individual for normalization into 'relative' load

varcounts <- data.frame()

for (sam in samples$sample) {
  high <- sum(df[df$impact == "HIGH", sam], na.rm = TRUE)
  high_hom <- sum(df[df$impact == "HIGH", sam] == 2, na.rm = TRUE) * 2
  mod <- sum(df[df$impact == "MODERATE", sam], na.rm = TRUE)
  mod_hom <- sum(df[df$impact == "MODERATE", sam] == 2, na.rm = TRUE) * 2
  low <- sum(df[df$impact == "LOW", sam], na.rm = TRUE)
  low_hom <- sum(df[df$impact == "LOW", sam] == 2, na.rm = TRUE) * 2
  alts <- sum(df[, sam], na.rm = TRUE)
  nsites <- sum(!is.na(df[, sam]), na.rm = TRUE)
  row <- c(sam, high, high_hom, mod, mod_hom, low, low_hom, alts, nsites)
  varcounts <- rbind(varcounts, row)
}

colnames(varcounts) <- c(
  "sample", "high", "high_hom", "mod", "mod_hom", "low", "low_hom", "alts",
  "nsites"
)
numcols <- c(
  "high", "high_hom", "mod", "mod_hom", "low", "low_hom", "alts",
  "nsites"
)
varcounts[numcols] <- sapply(varcounts[numcols], as.numeric)
varcounts <- merge(varcounts, samples, by = "sample")

write.table(varcounts,
  file = snakemake@output[["varcounts_nomiss"]], quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)