##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(meta.arrR)

source("05_Various/setup.R")

#### setup stuff ####

# number of local metaecosystems
n <- 5

# set min_per_i
min_per_i <- 120

# run the model for n years
years <- 5

max_i <- (60 * 24 * 365 * years) / min_per_i

# setup nutrient input to be maximum 10% of inital nutrients
input_mn <- meta.arrR::meta.arrR_starting_values$nutrients_pool * 0.1

#### Comparison  for increasing variability ####

# variability
variability <- seq(from = 0, to = 1, by = 0.25)

result <- purrr::map_dfr(variability, function(i) {
  
  meta.arrR::sim_nutr_input(n = n, max_i = max_i, variability = c(i, 0.0),
                            input_mn = input_mn, freq_mn = years) %>% 
    get_input_df() %>% 
    dplyr::mutate(Variability = i)}) %>% 
  tidyr::pivot_longer(-c(Timestep, Variability), names_to = "Meta", values_to = "Value") %>% 
  dplyr::filter(Meta != "Gamma") %>% 
  dplyr::mutate(Variability = paste0("Variability_", Variability))

#### Create ggplot ####

ggplot(data = result) + 
  geom_hline(yintercept = input_mn, linetype = 2, col = "grey") +
  geom_line(aes(x = Timestep, y = Value, col = Meta)) + 
  scale_color_viridis_d(name = "") + 
  facet_wrap(. ~ Variability, nrow = 1) + 
  theme_classic() + 
  theme(legend.position = "bottom")


##---------------------------##
#### Look at single result ####
##---------------------------##

# get data from HPC and save to ~/Downloads/results
filenames <- list.files(path = "02_Data/nofish/",
                        full.names = TRUE, pattern = "^meta-rn_nofish_*") %>%
  stringr::str_sort(numeric = TRUE)

# get df with variabilitities
variability_exp <- readRDS("02_Data/variability_exp.rds")

# add file names to df
variability_exp <- dplyr::bind_cols(variability_exp, file = filenames)

result_sub <- dplyr::slice(variability_exp, seq(from = 1, to = 200, by = 50)) %>% 
  dplyr::pull(file) %>% 
  purrr::map(readRDS)
