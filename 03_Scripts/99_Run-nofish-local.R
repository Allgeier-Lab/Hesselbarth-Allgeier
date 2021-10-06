##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

source("05_Various/setup.R")

#### Load data ####

default_starting <- readRDS("02_Data/default_starting.rds")

default_parameters <- readRDS("02_Data/default_parameters.rds")

#### Basic parameters ####

# set min_per_i
min_per_i <- 120

# run the model for n years
years <- 25

max_i <- (60 * 24 * 365 * years) / min_per_i

# run seagrass only 1 day
days <- 1

seagrass_each <- (24 / (min_per_i / 60)) * days

# save results only every m days
days <- 125 # which(max_i %% ((24 / (min_per_i / 60)) * (1:365)) == 0)

save_each <- (24 / (min_per_i / 60)) * days

max_i %% save_each

# set frequency of input peaks
freq_mn <- years * 1/4

# number of local metaecosystems
n <- 9

# create no reefs
reefs <- NULL

# setup extent and grain
dimensions <- c(100, 100)

grain <- 1

#### Adapt parameters ####

default_starting$pop_n <- 0

default_parameters$nutrients_diffusion <- 0.0

default_parameters$detritus_diffusion <- 0.0

default_parameters$detritus_fish_diffusion <- 0.0

#### Stable values ####

stable_values <- arrR::get_stable_values(starting_values = default_starting,
                                         parameters = default_parameters)

default_starting$nutrients_pool <- stable_values$nutrients_pool

default_starting$detritus_pool <- stable_values$detritus_pool

input_mn <- stable_values$nutr_input

#### Setup experiment ####

# # define parts to sample CV
# part <- c("bg_biomass", "ag_biomass", "bg_production", "ag_production")
# 
# # set low and high variability treatment levels
# variability_lo <- 0.0
# 
# variability_hi <- 0.5

#### Setup metaecosystems ####

# setup metaecosystems
metasyst <- meta.arrR::setup_meta(n = n, max_i = max_i, 
                                  starting_values = default_starting, parameters = default_parameters,
                                  dimensions = dimensions, grain = grain, reefs = reefs)

# simulate input 
input_temp <- meta.arrR::sim_nutr_input(n = n, max_i = max_i,
                                        variability = c(0.0, 0.0),
                                        input_mn = input_mn, freq_mn = freq_mn)

get_input_df(input_temp, gamma = FALSE) %>%
  dplyr::select(-Timestep) %>%
  apply(MARGIN = 2, FUN = sum)

stable_values$nutr_input * max_i

plot(input_temp, gamma = FALSE) + 
  geom_hline(yintercept = input_mn, linetype = 2, col = "black")

#### Run model ####

# run model
result <- meta.arrR::run_meta(metasyst = metasyst, nutr_input = input_temp, 
                              parameters = default_parameters,
                              max_i = max_i, min_per_i = min_per_i, 
                              seagrass_each = seagrass_each, save_each = save_each)

#### plot result ####

plot(result, summarize = TRUE)

plot_meta_production(result, lag = TRUE)

#### Calc variability ####

calc_variability(result, what = "biomass")

sample_variability(result, what = "biomass", lag = TRUE, verbose = FALSE)
