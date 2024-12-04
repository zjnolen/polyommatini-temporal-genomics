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

# Make it so we only use sites with 80% completeness in both modern and
# historical samples
comp80both <- (
  rowSums(!is.na(histgt)) >= floor(ncol(histgt) * 0.8) &
    rowSums(!is.na(modgt)) >= floor(ncol(modgt) * 0.8)
)

histgt <- histgt[comp80both, ]
modgt <- modgt[comp80both, ]
impacts <- impacts[comp80both]
df <- df[comp80both, ]

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
  inter <- sum(df[df$consequence == "intergenic_variant", sam], na.rm = TRUE)
  inter_hom <- sum(
    df[df$consequence == "intergenic_variant", sam] == 2,
    na.rm = TRUE
  ) * 2
  alts <- sum(df[, sam], na.rm = TRUE)
  nsites <- sum(!is.na(df[, sam]), na.rm = TRUE)
  row <- c(
    sam, high, high_hom, mod, mod_hom, low, low_hom, inter, inter_hom,
    alts, nsites
  )
  varcounts <- rbind(varcounts, row)
}

colnames(varcounts) <- c(
  "sample", "high", "high_hom", "mod", "mod_hom", "low", "low_hom", "inter",
  "inter_hom", "alts", "nsites"
)
numcols <- c(
  "high", "high_hom", "mod", "mod_hom", "low", "low_hom", "inter",
  "inter_hom", "alts", "nsites"
)
varcounts[numcols] <- sapply(varcounts[numcols], as.numeric)
varcounts <- merge(varcounts, samples, by = "sample")

write.table(varcounts,
  file = snakemake@output[["varcounts"]], quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)
