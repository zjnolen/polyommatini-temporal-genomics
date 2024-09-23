# taken from workflow/scripts/average_downsampled_thetas.R in
# https://zenodo.org/records/13383389

sink(file(snakemake@log[[1]], open = "wt"), type = "message")

library(Hmisc)

# Function that averages thetas across windows from sliding window output.
average_pestpg <- function(pestpg, popname, subsize, rep, minsites) {
  theta <- as.data.frame(read.table(pestpg, header = TRUE, comment.char = ""))
  theta <- theta[!theta$nSites < minsites, ]
  theta$watterson <- theta$tW / theta$nSites
  theta$pi <- theta$tP / theta$nSites
  watterson <- smean.cl.boot(
    theta[theta$nSites >= 1000, ]$watterson,
    B = 1000,
    na.rm = TRUE
  )
  pi <- smean.cl.boot(
    theta[theta$nSites >= 1000, ]$pi,
    B = 1000,
    na.rm = TRUE
  )
  tajima <- smean.cl.boot(
    theta[theta$nSites >= 1000, ]$Tajima,
    B = 1000,
    na.rm = TRUE
  )
  output <- theta[, c("Chr", "WinCenter", "pi", "watterson", "Tajima")]
  return(
    list(c(
      popname,
      subsize,
      rep,
      pi,
      watterson,
      tajima
    ),
    output)
  )
}

# averaging for input file
output <- average_pestpg(
  snakemake@input[[1]],
  snakemake@wildcards[["population"]],
  snakemake@wildcards[["samplesize"]],
  snakemake@wildcards[["rep"]],
  snakemake@params[["minsites"]]
)

# Write averages to file.
write.table(
  output[[1]],
  file = snakemake@output[[1]],
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  output[[2]],
  file = snakemake@output[[2]],
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)