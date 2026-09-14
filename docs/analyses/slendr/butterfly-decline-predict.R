library(slendr)
init_env()

args <- commandArgs(trailingOnly = TRUE)

# Script to make a slendr model based on the size changes inferred with gone.
# A time point in the past, present, and future can be sampled arbitrarily and
# the proportion retained heterozygosity relative to the past sampling point
# will be calculated for the past and future point. Additionally, runs of
# homozygosity will be estimated for each sample point, and the proportion of
# heterozygosity retained will be adjusted by these to assess how much RoH over
# a certain size contribute to the change. Despite being called 'butterfly',
# the basics of this script can be used with any gone trajectories, and despite
# the 'decline' they don't have to be declines. But it is kind of focused on
# seeing genetic diversity loss over a decline and seeing how it progresses if
# populations don't rebound. As of now the seq_len, rec_rate, and mut_rate are
# set to reasonable estimates for Lycaenid butterflies, so change those if you
# want to use it for something else. The minimum RoH size of 200 kb is also set
# to reflect ~ 200 generations of coalescence on average.

# Example usage:
# Rscript --vanilla butterfly-decline.R <gone-file> <mod-year> <hist-year> <fix-first-four> <max-change> <gen-time> <future-years>
# where:
# <gone-file> is the path to the gone file to base the model on
# <mod-year> is the year of the modern sample
# <hist-year> is the year of the historical sample
# <fix-first-four> is a boolean indicating whether the first four generations should be used as in the gone file or fixed to the value of the 5th
# <max-change> is the number of generations back in time to base the model off, the size at this generation will be used as the ancestral size
# <gen-time> is how many generations per year your organism has
# <future-years> how many years after <mod-year> to sample the 'future' values. The effective population sizes will be fixed to the <mod-year> sizes over these years.

gone <- as.character(args[1])
mod_year <- as.integer(args[2])
hist_year <- as.integer(args[3])
fix_first_four <- as.logical(args[4])
max_change <- as.integer(args[5])
gen_time <- as.integer(args[6])
extra_time <- as.integer(args[7])

if (fix_first_four) {
  gone <- read.table(gone, skip = 1, header = TRUE)
  gone[1:4, 2] <- gone[5, 2]
} else {
  gone <- read.table(gone, skip = 1, header = TRUE)
}

pop <- population("population", time = max_change + extra_time, N = gone[max_change, 2])
for (i in (max_change - 1):1) {
  pop <- resize(pop, time = i + extra_time, N = gone[i, 2], how = "step")
}

model <- compile_model(pop, generation_time = gen_time, direction = "backward")
samples <- schedule_sampling(model, times = c(mod_year-hist_year+extra_time, extra_time, 0), list(pop, 30))
seq_len <- 2e7
rec_rate <- 1.25e-8
mut_rate <- 2.9e-9
min_roh <- 200000

hz_froh_model_run <- function(model, seq_len, rec_rate, mut_rate, samples) {
  # Run the model
  ts <- msprime(model, sequence_length = seq_len, recombination_rate = rec_rate,
                          samples = samples, verbose = FALSE)
  # Overlay mutations
  ts <- ts_mutate(ts, mutation_rate = mut_rate)
  # Get sample list
  s <- ts_samples(ts)
  # grab samples by time period
  hist <- s[s$time == mod_year-hist_year+extra_time, ]$name
  mod <- s[s$time == extra_time, ]$name
  fut <- s[s$time == 0, ]$name
  # calculate the mean heterozygosity for each time period and the ratio
  hist_hz <- mean(ts_diversity(ts, hist)$diversity)
  mod_hz <- mean(ts_diversity(ts, mod)$diversity)
  fut_hz <- mean(ts_diversity(ts, fut)$diversity)
  mod_hist_hz_ratio <- mod_hz / hist_hz
  fut_hist_hz_ratio <- fut_hz / hist_hz
  # Identify runs of homozygosity in the samples
  ibd <- ts_ibd(ts, minimum_length = min_roh)
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
  froh_fut <- c()
  for (i in fut) {
    if (nrow(ibd[ibd$name1 == i & ibd$name2 == i, ]) == 0) {
      froh_fut <- c(froh_fut, 0)
    } else {
      froh_fut <- c(froh_fut, ibd[ibd$name1 == i & ibd$name2 == i, ]$total / seq_len)
    }
  }
  froh_fut <- mean(froh_fut)
  # Calculate heterozygosity adjusted for RoH (i.e. outside of RoH)
  hist_adj_hz <- mean(hist_hz / (1 - froh_hist))
  mod_adj_hz <- mean(mod_hz / (1 - froh_mod))
  fut_adj_hz <- mean(fut_hz / (1 - froh_fut))
  mod_hist_adj_hz_ratio <- mod_adj_hz / hist_adj_hz
  fut_hist_adj_hz_ratio <- fut_adj_hz / hist_adj_hz
  # Output results
  return(round(c(mod_hist_hz_ratio, fut_hist_hz_ratio, mod_hist_adj_hz_ratio, fut_hist_adj_hz_ratio, froh_hist, froh_mod, froh_fut), 5))
}

hz_froh_model_run(model, seq_len, rec_rate, mut_rate, samples)
