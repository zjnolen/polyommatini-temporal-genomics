# taken from workflow/scripts/average_downsampled_thetas.R in
# https://zenodo.org/records/13383389

sink(file(snakemake@log[[1]], open = "wt"), type = "message")

# Function that averages thetas across windows from sliding window output.
average_pestpg <- function(pestpg, popname, subsize, rep, minsites) {
  theta <- as.data.frame(read.table(pestpg, header = TRUE, comment.char = ""))
  theta <- theta[!theta$nSites < minsites, ]
  theta$watterson <- theta$tW / theta$nSites
  theta$pi <- theta$tP / theta$nSites
  return(
    cbind(
      popname,
      subsize,
      rep,
      mean(theta$pi),
      mean(theta$watterson),
      mean(theta$Tajima)
    )
  )
}

# averaging for input files
avg <- average_pestpg(
  snakemake@input[[1]],
  snakemake@wildcards[["population"]],
  snakemake@wildcards[["samplesize"]],
  snakemake@wildcards[["rep"]],
  snakemake@params[["minsites"]]
)

# Write averages to file.
write.table(
  avg,
  file = snakemake@output[[1]],
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)