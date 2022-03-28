##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

source("05_Various/setup.R")

#### Load data ####

starting_list <- readRDS("02_Data/starting_list.rds")

parameters_list <- readRDS("02_Data/parameters_list.rds")

variability_experiment <- readRDS("02_Data/variability_experiment.rds")

variability_experiment[c(1, 51, 101, 151, 201, 251, 301, 351, 401), ]

meta_rn <- readRDS("02_Data/nofish/meta-rn_nofish_201.rds")

#### Plot input and output ####

plot(meta_rn$nutr_input, gamma = FALSE)

plot(meta_rn, summarize = TRUE)

#### Plot production #### 

meta_production <- purrr::map_dfr(meta_rn$seafloor, function(i) {
  
  dplyr::select(i, x, y, timestep, bg_production, ag_production) %>% 
    dplyr::group_by(timestep) %>% 
    dplyr::summarise(bg_prod_sum = sum(bg_production), 
                     ag_prod_sum = sum(ag_production), .groups = "drop") %>% 
    dplyr::mutate(bg_diff = bg_prod_sum - dplyr::lag(bg_prod_sum), 
                  ag_diff = ag_prod_sum - dplyr::lag(ag_prod_sum))}, .id = "meta") %>% 
  tidyr::pivot_longer(-c(meta, timestep))

dplyr::filter(meta_production, name %in% c("bg_prod_sum", "ag_prod_sum")) %>% 
  ggplot() + 
  geom_line(aes(x = timestep, y = value, col = meta)) + 
  facet_wrap(. ~ name, scales = "free_y", ncol = 1) + 
  scale_color_viridis_d(option = "C") + 
  theme_classic() + 
  theme(legend.position = "none")

dplyr::filter(meta_production, name %in% c("bg_diff", "ag_diff")) %>% 
  ggplot() + 
  geom_line(aes(x = timestep, y = value, col = meta)) + 
  facet_wrap(. ~ name, scales = "free_y", ncol = 1) + 
  scale_color_viridis_d(option = "C") + 
  theme_classic() + 
  theme(legend.position = "none")
