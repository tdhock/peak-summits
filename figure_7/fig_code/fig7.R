options(namedCapture.engine="PCRE")
works_with_R(
  "3.5.1",
  penaltyLearning="2018.9.4",
  data.table="1.11.6",
  ggplot2="3.0.0",
  tidyverse="1.2.1",
  latex2exp="0.4.0",
  "jewellsean/FastLZeroSpikeInference@e4fcf32333bac8d01122b5352032b1aa28c6f714")

source("utils.R")

## data parameters
fps <- 100 # frames per second 
data_set <- 7 # spike finder dataset number 
ind <- 13 # cell number 
subset_is <- 1:20000 # subset of the data used to estimate the solution 
gam <- 0.98 # decay parameter 

## tuning parameters 
lambda_p2 <- 3.6
lambda_p3 <- 3.76

## file i/o paramters 
spikefinder_base_dir <- "../spike_finder_data/"
fig_save_dir <- "../figures/"

## 1. Read and process spike finder data
data_dir <- paste0(spikefinder_base_dir, "spikefinder.train/", data_set, ".train.")
calcium_dat <- readr::read_csv(paste0(data_dir, "calcium.csv")) %>% rename_all(col_renamer)
spike_dat <- readr::read_csv(paste0(data_dir, "spikes.csv")) %>% rename_all(col_renamer)
data <- subset_data(calcium_dat, spike_dat, subset_is)

## 2. Estimate spikes and calcium concentrations corresponding to eqn. (2) and (3)
fit_2 <- estimate_spikes(dat = data$calcium_d, gam = gam, lambda = lambda_p2, 
                         constraint = T, estimate_calcium = T)


times <- subset_is * (1 / fps)
dt <- data.table(seconds=times, calcium=data$calcium_d, AR1=fit_2$estimated_calcium)
fwrite(dt[, list(calcium)], "../../Cpattern1D/data.txt", col.names=FALSE)

if(FALSE){
  neuroSpike <- dt
  save(neuroSpike, file="~/R/gfpop/data/neuroSpike.RData", compress="xz")
  prompt(neuroSpike, file="~/R/gfpop/man/neuroSpike.Rd")
}

