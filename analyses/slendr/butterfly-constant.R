library(slendr)
init_env()

args <- commandArgs(trailingOnly = TRUE)

# Script to make a slendr model for sampling two time points under an assumed
# constant population size. This gives an idea of the expected heterozygosity
# and Froh in a population at a given equilibrium size and also the amount of
# fluctuations expectable in these values over a sample period. Despite the
# name, it is only butterfly oriented because of the current rec_rate and
# mut_rate settings, so can be used for other systems if you change those. Also
# assumes 1 generation = 1 year, so if that's not true for what you want to do
# just set mod-year to 0 and hist-year to the number of generations in the past
# you want to sample from.

# Example usage:
# Rscript --vanilla butterfly-constant.R <mod-year> <hist-year> <size> <max-time>
# where:
# <mod-year> is the year of the modern sample
# <hist-year> is the year of the historical sample
# <size> is the constant population size to simulate
# <max-time> is the number of generations back to mark the ancestral. This is
#  only really important in that it should be longer than the time between your
#  values for mod-year and hist-year so that the sampling points happen inside
#  the model definition. You can set it to a really big number as well, the
#  model will simulate backwards in time anyway to coalescense so it won't
#  affect the results.
#  
#  Will print output to standard out, be sure to direct it somewhere. Can run
#  multiple reps pretty quickly with GNU Parallel, e.g.:
#  parallel --jobs <threads> -N0 Rscript --vanilla butterfly-constant.R <opts> ::: {1..100} > model_output.txt



mod_year <- as.integer(args[1])
hist_year <- as.integer(args[2])
size <- as.integer(args[3])
max_time <- as.integer(args[4])

pop <- population("population", time = max_time, N = size)

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
