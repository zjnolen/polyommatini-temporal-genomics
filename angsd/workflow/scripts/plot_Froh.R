sink(file(snakemake@log[[1]], open = "wt"), type = "message")

library(dplyr)
library(ggplot2)

# Calculate Froh for separate user defined size bins. Describes the proportion
# of total Froh that comes from each size class
froh_bins <- function(samplelist, roh_df, minroh, bins, lenautos) {
  samples <- as.data.frame(read.table(samplelist, header = TRUE, sep = "\t"))

  norun <- data.frame(matrix(nrow = nrow(samples), ncol = 0))
  norun$group <- "RG"
  norun$sample <- samples$sample
  norun$chr <- rep(0, nrow(samples))
  norun$start <- rep(0, nrow(samples))
  norun$end <- rep(0, nrow(samples))
  norun$length <- rep(0, nrow(samples))
  norun$nsnps <- rep(0, nrow(samples))
  norun$phred <- rep(100, nrow(samples))

  bins <- c(minroh, bins, Inf)
  frohs <- c()

  for (i in seq_along(bins)) {
    if (i < length(bins)) {
      runs <- roh_df[roh_df$length >= bins[i], ]
      runs <- runs[runs$length < bins[i + 1], ]
      runs <- rbind(runs, norun)
      runs <- runs %>%
        group_by(sample) %>%
        summarize(length = sum(length))
      runs <- merge(runs, samples, by = "sample")
      runs$froh <- runs$length / as.numeric(lenautos)
      runs$range <- paste0("[", bins[i], ",", bins[i + 1], ")")
      frohs <- rbind(frohs, runs)
    }
  }

  return(frohs)
}

# Calculate the number of rohs and the cumulative length per individual for
# rohs > the minimum roh size. This can help distinguish patterns of
# consanguinity from increased background relatedness
nroh_cumroh <- function(roh_df, minroh, samplelist) {
  samples <- as.data.frame(read.table(samplelist, header = TRUE, sep = "\t"))
  df <- roh_df[roh_df$length >= minroh, ]
  df <- df %>%
    group_by(sample) %>%
    summarize(
      cumlen = sum(length),
      nroh = n()
    )
  df <- merge(df, samples, by = "sample")
  return(df)
}

# Plot binned Froh per individual
plot_bins <- function(frohs, minroh, bins, plotpre) {
  bins <- c(minroh, bins, Inf)
  ranges <- c()
  for (i in seq_along(bins)) {
    if (i < length(bins)) {
      ranges <- c(ranges, paste0("[", bins[i], ",", bins[i + 1], ")"))
    }
  }
  frohs$range <- as.factor(frohs$range)
  frohs$range <- factor(frohs$range, levels = rev(ranges), labels = rev(ranges))
  indfroh <- frohs %>%
    group_by(sample) %>%
    summarize(froh = sum(froh))

  frohs$sample <- factor(
    frohs$sample,
    levels = indfroh$sample[order(indfroh$froh)]
  )

  ggplot(data = frohs, aes(x = sample, y = froh, fill = range)) +
    geom_bar(position = "stack", stat = "identity") +
    facet_grid(cols = vars(population), scales = "free_x") +
    theme_bw() +
    labs(y = expression(FRoH), fill = "RoH Length (bp)") +
    scale_fill_grey() +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    theme(
      axis.title.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  ggsave(paste0(plotpre, ".froh_bins.svg"), width = 18, height = 6)
}

# Plot Nroh ~ Cumulative roh length for individuals with rohs
plot_cumroh_nroh <- function(nroh_cumroh, plotpre) {
  ggplot(data = nroh_cumroh, aes(x = cumlen, y = nroh, color = population)) +
    geom_point() +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    ylab("Total number of runs of homozygosity") +
    xlab("Total length of runs of homozygosity (bp)") +
    theme_classic() +
    labs(color = "Population")

  ggsave(paste0(plotpre, ".cumroh_nroh.svg"), width = 9, height = 7)
}


# Actually run functions to produce outputs, saving tables for the bed file of
# all rohs for all individuals and the Froh calculated for all rohs > the
# minimum roh size for each individual

if (file.size(snakemake@input[["roh"]]) > 0) {
  df <- read.table(snakemake@input[["roh"]], header = FALSE)
} else {
  df <- data.frame(matrix(ncol = 8, nrow = 0))
}

colnames(df) <- c(
  "group", "sample", "chr", "start", "end", "length", "nsnps", "phred"
)
df <- df[df$phred >= snakemake@params[["roh_phred"]], ]


frohs <- froh_bins(
  snakemake@input[["inds"]],
  df,
  snakemake@params[["minroh"]],
  snakemake@params[["bins"]],
  read.table(snakemake@input[["autos"]], sep = "\t")[, 2]
)

indfroh <- frohs %>%
  group_by(sample) %>%
  summarize(Froh = sum(froh))

indfroh <- merge(
  as.data.frame(
    read.table(snakemake@input[["inds"]], header = TRUE, sep = "\t")
  ),
  indfroh,
  by = "sample"
)

write.table(indfroh,
  file = paste0(snakemake@params[["outpre"]], ".ind_froh.tsv"), quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = TRUE
)

cumroh <- nroh_cumroh(
  df,
  snakemake@params[["minroh"]],
  snakemake@input[["inds"]]
)

plot_bins(
  frohs,
  snakemake@params[["minroh"]],
  snakemake@params[["bins"]],
  snakemake@params[["outpre"]]
)

plot_cumroh_nroh(
  cumroh,
  snakemake@params[["outpre"]]
)
