##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

#### Load setupt ####

source("05_Various/setup.R")

#### Adapt parameters ####

# parameters_list$seagrass_thres <- -1/2

parameters_list$nutrients_diffusion <- 0.0
parameters_list$detritus_diffusion <- 0.0
parameters_list$detritus_fish_diffusion <- 0.0

n <- 5

starting_list$pop_n <- 0

starting_list$bg_biomass <- parameters_list$bg_biomass_max

starting_list$ag_biomass <- parameters_list$ag_biomass_max

#### Stable values ####

stable_values <- arrR::get_stable_values(bg_biomass = starting_list$bg_biomass,
                                         ag_biomass = starting_list$ag_biomass,
                                         parameters = parameters_list)

starting_list$nutrients_pool <- stable_values$nutrients_pool

starting_list$detritus_pool <- stable_values$detritus_pool

input_mn <- stable_values$nutr_input

#### Setup metaecosystems ####

# setup metaecosystems
metasyst <- meta.arrR::setup_meta(n = n, max_i = max_i, 
                                  starting_values = starting_list, parameters = parameters_list,
                                  dimensions = dimensions, grain = grain, reef = NULL)

# simulate input 
input_temp <- meta.arrR::sim_nutr_input(n = n, max_i = max_i,
                                        # variability = c(1.0, 0.0),
                                        amplitude_mod = c(0, 1/4, 1/2, 3/4, 1.0),
                                        input_mn = input_mn, freq_mn = freq_mn)

# get_input_df(input_temp, gamma = FALSE) %>%
#   dplyr::select(-Timestep) %>%
#   apply(MARGIN = 2, FUN = sum)
# 
# stable_values$nutr_input * max_i

# input_temp$values <- purrr::map2(input_temp$values, input_temp$amplitude_i,
#                                  function(x,y) x + y)

# plot(input_temp, gamma = FALSE) +
#   geom_hline(yintercept = input_mn, linetype = 2, col = "black") +
#   geom_hline(yintercept = stable_values$nutr_input, linetype = 2, col = "grey")

#### Run model ####

# run model
result <- meta.arrR::run_meta(metasyst = metasyst, nutr_input = input_temp, 
                              parameters = parameters_list,
                              max_i = max_i, min_per_i = min_per_i, 
                              seagrass_each = seagrass_each, save_each = save_each)

#### plot result ####

plot(result, summarize = TRUE)

plot_meta_production(result, lag = TRUE) + 
  geom_vline(xintercept = 17700, linetype = 2) + 
  geom_vline(xintercept = 35100, linetype = 2)

#### Random stuff ####

# prod <- get_meta_production(result, lag = TRUE)
# 
# dplyr::filter(prod, timestep > 100000, timestep <= 200000) %>% 
#   dplyr::group_by(meta, part) %>% 
#   dplyr::summarise(value = sum(value), .groups = "drop") %>% 
#   dplyr::group_by(part) %>% 
#   dplyr::mutate(value = (value - max(value)) / max(value) * 100) %>% 
#   dplyr::arrange(part, meta)

purrr::map_dfr(result$seafloor, function(i) {
  
  dplyr::filter(i, timestep == max(timestep)) %>% 
    dplyr::summarise(ag_prod = sum(ag_production), 
                     bg_prod = sum(bg_production))}) %>% 
  cbind(amp = c(0, 1/4, 1/2, 3/4, 1.0), .) %>% 
  dplyr::mutate(ttl_prod = ag_prod + bg_prod, 
                ag_rel = (ag_prod - max(ag_prod)) / max(ag_prod) * 100, 
                bg_rel = (bg_prod - max(bg_prod)) / max(bg_prod) * 100, 
                ttl_rel = (ttl_prod - max(ttl_prod)) / max(ttl_prod) * 100)

# #### Calc variability ####
# variability <- sample_variability(x = result, what = "production", itr = 100)
# 
# dplyr::filter(variability, stat %in% c("gamma", "synchrony")) %>% 
#   ggplot() + 
#   geom_line(aes(x = factor(n), y = mean, group = stat), col = "grey") + 
#   geom_point(aes(x = factor(n), y = mean, col = stat)) + 
#   geom_linerange(aes(x = factor(n), ymin = mean - sd, ymax = mean + sd, col = stat)) + 
#   facet_wrap(. ~ part, nrow = 2) +
#   scale_y_continuous(limits = c(0, 1)) + 
#   theme_classic()
