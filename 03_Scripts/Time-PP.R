##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: 

#### Load setup ####

source("05_Various/setup.R")
source("01_Functions/import_data.R")

#### Get data ####

n <- 5

experiment_df <- readRDS("02_Data/00_experiment_df.rds")

noise_df <- readRDS(file = paste0("02_Data/05-variability-noise-", n, ".rds"))

#### Classifiy variability ####

experiment_df <- dplyr::mutate(experiment_df, biotic_class = cut(biotic, breaks = 3, labels = c("low", "med", "hi")), 
                                abiotic_class = cut(abiotic, breaks = 3, labels = c("low", "med", "hi"))) %>% 
  tidyr::unite("combined_class", c(biotic_class, abiotic_class))

experiment_class_df <- dplyr::group_by(experiment_df, combined_class) %>% 
  dplyr::group_split() %>% 
  purrr::map_dfr(function(x) {
    
    dplyr::mutate(x, biotic_mn = mean(biotic), abiotic_mn = mean(abiotic), 
                  biotic_diff = abs(biotic_mn - biotic), abiotic_diff =  abs(abiotic_mn - abiotic), 
                  total_diff = biotic_diff + abiotic_diff) %>% 
      dplyr::arrange(total_diff) %>% 
      head(n = 1)}) %>% 
  dplyr::arrange(biotic, dplyr::desc(abiotic)) %>% 
  dplyr::mutate(label = LETTERS[1:dplyr::n()])

sample_id <- which(experiment_df$biotic %in% experiment_class_df$biotic & 
                     experiment_df$abiotic %in% experiment_class_df$abiotic)

size_base <- 10

ggplot() + 
  geom_point(data = experiment_df, aes(x = biotic, y = abiotic)) + 
  geom_label(data = experiment_class_df, aes(x = biotic, y = abiotic, label = combined_class)) + 
  coord_equal() + 
  scale_color_viridis_d() +
  theme_classic(base_size = size_base) + 
  theme(legend.position = "none")

#### ####

noise_df[sample_id] %>% 
  purrr::set_names(experiment_df[sample_id, "combined_class", drop = TRUE]) %>% 
  purrr::map_dfr(function(x) x$prod_time, .id = "combined_class") %>% 
  dplyr::group_by(combined_class, timestep, pop_n) %>% 
  dplyr::summarise(ttl_production = sum(ttl_production)) %>% 
  ggplot(aes(x = timestep, y = ttl_production, color = factor(pop_n))) + 
  geom_line() +
  facet_wrap(. ~ combined_class, ncol = 3, nrow = 3) + 
  scale_color_viridis_d() + 
  theme_classic() + 
  theme(legend.position = "bottom")
