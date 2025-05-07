library(slendr)
init_env()

args <- commandArgs(trailingOnly = TRUE)

# Script to make a slendr model based on the size changes inferred with gone,
# then sample at two time points, calculating the proportion heterozygosity
# retained, the proportion of heterozygosity retained outside of runs of
# homozygosity, and the mean Froh for historical and modern samples. Despite
# the name, it can work on any gone output, not just butterflies, and would
# model whatever changes are in the gone output, not just declines. Just
# change the rec_rate and mut_rate, which are typical of butterflies for now.
# Also assumes 1 generation = 1 year, so if that's not true for what you want to
# do just set mod-year to 0 and hist-year to the number of generations in the
# past you want to sample from.

# Example usage:
# Rscript --vanilla butterfly-decline.R <gone-file> <mod-year> <hist-year> <fix-first-four> <max-change>
# where:
# <gone-file> is the path to the gone file to base the model on
# <mod-year> is the year of the modern sample
# <hist-year> is the year of the historical sample
# <fix-first-four> is a boolean indicating whether the first four generations should be used as in the gone file or fixed to the value of the 5th
# <max-change> is the number of generations back in time to base the model off, the size at this generation will be used as the ancestral size
#  
#  Will print output to standard out, be sure to direct it somewhere. Can run
#  multiple reps pretty quickly with GNU Parallel, e.g.:
#  parallel --jobs <threads> -N0 Rscript --vanilla butterfly-decline.R <opts> ::: {1..100} > model_output.txt

gone <- as.character(args[1])
mod_year <- as.integer(args[2])
hist_year <- as.integer(args[3])
fix_first_four <- as.logical(args[4])
max_change <- as.integer(args[5])

if (fix_first_four) {
  gone <- read.table(gone, skip = 1, header = TRUE)
  gone[1:4, 2] <- gone[5, 2]
} else {
  gone <- read.table(gone, skip = 1, header = TRUE)
}

pop <- population("population", time = max_change, N = gone[max_change, 2])
for (i in (max_change - 1):1) {
  pop <- resize(pop, time = i, N = gone[i, 2], how = "step")
}

model <- compile_model(pop, generation_time = 1, direction = "backward")
samples <- schedule_sampling(model, times = c(mod_year - hist_year, 0), list(pop, 30))
seq_len <- 1e7
rec_rate <- 5e-9
mut_rate <- 2.9e-8

hz_froh_model_run <- function(model, seq_len, rec_rate, mut_rate, samples) {
  # Run the model
  ts <- msprime(model,
    sequence_length = seq_len, recombination_rate = rec_rate,
    samples = samples, verbose = FALSE
  )
  # Overlay mutations
  ts <- ts_mutate(ts, mutation_rate = mut_rate)
  # Get sample list
  s <- ts_samples(ts)
  # grab samples by time period
  hist <- s[s$time == mod_year - hist_year, ]$name
  mod <- s[s$time == 0, ]$name
  # calculate the mean heterozygosity for each time period and the ratio
  hist_hz <- mean(ts_diversity(ts, hist)$diversity)
  mod_hz <- mean(ts_diversity(ts, mod)$diversity)
  mod_hist_hz_ratio <- mod_hz / hist_hz
  # Identify runs of homozygosity in the samples
  ibd <- ts_ibd(ts, minimum_length = 100000)
  # Calculate the mean IBD for each time period
  froh_hist <- c()
  for (i in hist) {
    if (nrow(ibd[ibd$name1 == i & ibd$name2 == i, ]) == 0) {
      froh_hist <- c(froh_hist, 0)
    } else {
      froh_hist <- c(froh_hist, ibd[ibd$name1 == i & ibd$name2 == i, ]$total / seq_len)
    }
  }
  froh_hist <- mean(froh_hist)
  froh_mod <- c()
  for (i in mod) {
    if (nrow(ibd[ibd$name1 == i & ibd$name2 == i, ]) == 0) {
      froh_mod <- c(froh_mod, 0)
    } else {
      froh_mod <- c(froh_mod, ibd[ibd$name1 == i & ibd$name2 == i, ]$total / seq_len)
    }
  }
  froh_mod <- mean(froh_mod)
  # Calculate heterozygosity adjusted for RoH (i.e. outside of RoH)
  hist_adj_hz <- mean(hist_hz / (1 - froh_hist))
  mod_adj_hz <- mean(mod_hz / (1 - froh_mod))
  mod_hist_adj_hz_ratio <- mod_adj_hz / hist_adj_hz
  # Output results
  return(c(mod_hist_hz_ratio, mod_hist_adj_hz_ratio, froh_hist, froh_mod))
}

hz_froh_model_run(model, seq_len, rec_rate, mut_rate, samples)
