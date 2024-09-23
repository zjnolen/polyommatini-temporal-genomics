sink(file(snakemake@log[[1]], open = "wt"), type = "message")

library(rcompanion) # to get Z-score from test

nreps <- as.integer(snakemake@params[["nreps"]])
subsample_n <- as.integer(snakemake@wildcards[["samplesize"]])
pop1 <- snakemake@wildcards[["population1"]]
pop2 <- snakemake@wildcards[["population2"]]


stats <- c()

for (rep in c(1:nreps)) {
  p1 <- read.table(snakemake@input[["pop1"]][rep], header = TRUE, sep = "\t")
  p2 <- read.table(snakemake@input[["pop2"]][rep], header = TRUE, sep = "\t")
  df <- merge(p1, p2, by = c("Chr", "WinCenter"))
  nwins <- nrow(df)
  pidiff <- median(df$pi.x - df$pi.y)
  wattersondiff <- median(df$watterson.x - df$watterson.y)
  tajimadiff <- median(df$Tajima.x - df$Tajima.y)
  piret <- median(df$pi.x / df$pi.y)
  wattersonret <- median(df$watterson.x / df$watterson.y)
  tajimaret <- median(df$Tajima.x / df$Tajima.y)
  pval_pi <- wilcox.test(df$pi.x, df$pi.y, paired = TRUE)$p.value
  pval_watterson <- wilcox.test(
    df$watterson.x, df$watterson.y,
    paired = TRUE
  )$p.value
  pval_tajima <- wilcox.test(df$Tajima.x, df$Tajima.y, paired = TRUE)$p.value
  z_pi <- wilcoxonZ(df$pi.x, df$pi.y, paired = TRUE)[[1]]
  z_watterson <- wilcoxonZ(df$watterson.x, df$watterson.y, paired = TRUE)[[1]]
  z_tajima <- wilcoxonZ(df$Tajima.x, df$Tajima.y, paired = TRUE)[[1]]
  effsize_pi <- z_pi / sqrt(nwins)
  effsize_watterson <- z_watterson / sqrt(nwins)
  effsize_tajima <- z_tajima / sqrt(nwins)
  row <- c(
    pop1, pop2, subsample_n, rep,
    pidiff, piret, z_pi, effsize_pi, pval_pi,
    wattersondiff, wattersonret, z_watterson, effsize_watterson, pval_watterson,
    tajimadiff, tajimaret, z_tajima, effsize_tajima, pval_tajima
  )
  stats <- rbind(stats, row)
}

stats <- data.frame(stats)

write.table(
  stats,
  file = snakemake@output[["stats"]],
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)
