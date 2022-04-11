# source setup file
source("05_Various/setup.R")

source("01_Functions/update_starting.R")

nutr_req <- function(bg_biomass, ag_biomass, parameters) {
  
  # calculate detritus modifier for bg biomass
  bg_modf <- (bg_biomass - parameters$bg_biomass_min) /
    (parameters$bg_biomass_max - parameters$bg_biomass_min)
  
  # calculate detritus modifier for ag biomass
  ag_modf <- (ag_biomass - parameters$ag_biomass_min) /
    (parameters$ag_biomass_max - parameters$ag_biomass_min)
  
  # calculate ag detritus
  bg_detritus <- bg_biomass * (parameters$seagrass_slough * bg_modf)
  
  # calculate ag detritus
  ag_detritus <- ag_biomass * (parameters$seagrass_slough * ag_modf)
  
  # MH: nutrients pool neither includes sloughed nutrients nor output
  
  # combine to list
  list(bg = bg_detritus * parameters$bg_gamma, ag = ag_detritus * parameters$ag_gamma)
  
}

get_lag_slough <- function(result) {
  
  slough_temp <- lapply(X = seq_along(result$seafloor), FUN = function(i) {
    
    # select only required columns
    seafloor_temp <- subset(x = result$seafloor[[i]],
                            select = c("timestep", "ag_slough", "bg_slough"))
    
    # sum for each timestep
    seafloor_temp <- stats::aggregate(x = seafloor_temp[, -1],
                                      by = list(timestep = seafloor_temp$timestep),
                                      FUN = "sum")
    
    # ag_production
    seafloor_temp[, "ag_slough"] <-
      c(NA, seafloor_temp[2:nrow(seafloor_temp), "ag_slough"] -
          seafloor_temp[1:(nrow(seafloor_temp) - 1), "ag_slough"])
    
    # bg_production
    seafloor_temp[, "bg_slough"] <-
      c(NA, seafloor_temp[2:nrow(seafloor_temp), "bg_slough"] -
          seafloor_temp[1:(nrow(seafloor_temp) - 1), "bg_slough"])
      
    # reshape to long format
    seafloor_temp <- stats::reshape(data = seafloor_temp, direction = "long",
                                    v.names = "value", varying = list(names(seafloor_temp[, -1])),
                                    idvar = "timestep", ids = seafloor_temp[, 1],
                                    timevar = "part", times = names(seafloor_temp[, -1]),
                                    new.row.names = seq(from = 1, to = nrow(seafloor_temp) * 2))
    
    # timestep not needed here, will be added in cbind call
    cbind(meta = i, seafloor_temp)
    
  })
  
  # combine to one data.frame
  slough_temp <- do.call(what = "rbind", args = slough_temp)
  
  return(slough_temp)
  
}

# load parameters
starting_list <- readRDS("02_Data/starting_list.rds")

parameters_list <- readRDS("02_Data/parameters_list.rds")

#### Adapt parameters ####

starting_list$pop_n <- 0

parameters_list$nutrients_diffusion <- 0.0
parameters_list$detritus_diffusion <- 0.0
parameters_list$detritus_fish_diffusion <- 0.0

n <- 1

# parameters_list$seagrass_thres <- 1/3

#### Update starting ####

# update_starting(parameters = parameters_list, starting = starting_list, 
#                 ag = 1/3, bg = 0.95)

#### Stable values ####

stable_values <- arrR::get_stable_values(bg_biomass = starting_list$bg_biomass,
                                         ag_biomass = starting_list$ag_biomass,
                                         parameters = parameters_list)

starting_list$nutrients_pool <- stable_values$nutrients_pool

starting_list$detritus_pool <- stable_values$detritus_pool

#### Setup input and metaecosyst ####

# simulate input
input_temp <- meta.arrR::sim_nutr_input(n = n, max_i = max_i, freq_mn = freq_mn,
                                        input_mn = stable_values$nutr_input * 0.5,
                                        amplitude_mod = runif(n = n, min = 0.0, max = 0),
                                        phase_mod = runif(n = n, min = 0.0, max = 0))

# plot(input_temp, gamma = FALSE) +
#   geom_hline(yintercept = stable_values$nutr_input, linetype = 2)

# setup metaecosystems
metasyst_temp <- meta.arrR::setup_meta(n = n, max_i = max_i,
                                       starting_values = starting_list,
                                       parameters = parameters_list,
                                       dimensions = dimensions, grain = grain,
                                       reef = NULL)

# run model
result_temp <- meta.arrR::run_meta(metasyst = metasyst_temp, nutr_input = input_temp,
                                   parameters = parameters_list,
                                   max_i = max_i, min_per_i = min_per_i,
                                   seagrass_each = seagrass_each, save_each = save_each)

# result_temp <- readRDS("02_Data/result_temp.rds")
# input_temp <- result_temp$nutr_input

# result_temp <- meta.arrR::filter_meta(x = result_temp, filter = c(max_i / 2, max_i))

plot(result_temp, summarize = TRUE)
plot_meta_production(result_temp, lag = TRUE)

dplyr::filter(result_temp$seafloor$Meta_1, timestep == 0) %>% 
  dplyr::select(timestep, ag_biomass, bg_biomass, nutrients_pool, detritus_pool)





parameters_list$bg_biomass_min + ((parameters_list$bg_biomass_max - 
parameters_list$bg_biomass_min) * 0.5)

result_sum <- purrr::map_dfr(result_temp$seafloor, function(i) {
  
  dplyr::group_by(i, timestep) %>% 
    dplyr::summarise(bg_biomass = sum(bg_biomass), ag_biomass = sum(ag_biomass),
                     bg_production = sum(bg_production), ag_production = sum(ag_production), 
                     .groups = "drop") %>% 
    dplyr::mutate(ttl_biomass = bg_biomass + ag_biomass, 
                  ttl_production = bg_production + ag_production)}, .id = "meta")

result_sum_long <- dplyr::select(result_sum, meta, timestep, ag_biomass, bg_biomass) %>% 
  tidyr::pivot_longer(-c(meta, timestep), names_to = "part") %>%
  dplyr::arrange(part, meta, timestep)

result_prod <- meta.arrR::get_meta_production(result_temp, lag = TRUE)

result_slough <- get_lag_slough(result = result_temp)

#### Biomass and Slough amounts ###

# total nutrients pool
purrr::map_dfr(result_temp$seafloor, function(i) {
  
  dplyr::select(i, timestep, nutrients_pool, ag_uptake, bg_uptake) %>% 
    dplyr::group_by(timestep) %>% 
    dplyr::summarise(ag_uptake = sum(ag_uptake), bg_uptake = sum(bg_uptake), 
                     nutrients_pool = sum(nutrients_pool), .groups = "drop") %>% 
    dplyr::mutate(ttl_uptake = ag_uptake + bg_uptake) %>% 
    dplyr::select(-c(ag_uptake, bg_uptake)) %>%
    tidyr::pivot_longer(-timestep)}, .id = "meta") %>% 
  dplyr::filter(name == "nutrients_pool") %>% 
  ggplot() + 
  geom_line(aes(x = timestep, y = value, col = meta)) + 
  facet_wrap(. ~ name, scales = "free_y") + 
  scale_color_viridis_d(name = "", option = "C") +
  labs(x = "Timestep", y = "Total nutrients pool") + 
  guides(col = guide_legend(nrow = 1)) +
  theme_classic() + theme(legend.position = "bottom")

# total biomass values
dplyr::select(result_sum, meta, timestep, bg_biomass, ag_biomass) %>% 
  tidyr::pivot_longer(-c(meta, timestep)) %>% 
  dplyr::mutate(name = factor(name, levels = c("ag_biomass", "bg_biomass"))) %>% 
  ggplot() + 
  geom_line(aes(x = timestep, y = value, col = meta)) + 
  facet_wrap(. ~ name, scales = "free_y", ncol = 1) + 
  scale_color_viridis_d(name = "", option = "C") +
  labs(x = "Timestep", y = "Total biomass") + 
  guides(col = guide_legend(nrow = 1)) +
  theme_classic() + theme(legend.position = "bottom")

# difference total production slough
purrr::map_dfr(result_temp$seafloor, function(i) {
  
  dplyr::select(i, timestep, ag_production, bg_production, 
                ag_slough, bg_slough) %>% 
    dplyr::group_by(timestep) %>% 
    dplyr::summarise(ag_prod = sum(ag_production), bg_prod = sum(bg_production), 
                     ag_slo = sum(ag_slough), bg_slo = sum(bg_slough), .groups = "drop") %>% 
    dplyr::mutate(ag_difference = ag_prod - ag_slo, 
                  bg_difference = bg_prod - bg_slo) %>% 
    dplyr::select(timestep, ag_difference, bg_difference) %>% 
    tidyr::pivot_longer(-timestep) %>% 
    dplyr::mutate(name = factor(name, levels = c("ag_difference", "bg_difference")))}, .id = "meta") %>% 
  ggplot() + 
  geom_line(aes(x = timestep, y = value, col = meta)) + 
  facet_wrap(. ~ name, scales = "free_y", ncol = 1) + 
  # geom_hline(yintercept = 0, linetype = 2) + 
  scale_color_viridis_d(name = "", option = "C") +
  labs(x = "Timestep", y = "Total production - Total slough") + 
  guides(col = guide_legend(nrow = 1)) +
  theme_classic() + theme(legend.position = "bottom")

# difference cumulative production slough
dplyr::bind_cols(result_prod, 
                 part.slough = result_slough$part, value.slough = result_slough$value, 
                 part.biom = result_sum_long$part, value.biom = result_sum_long$value) %>% 
  dplyr::mutate(value.diff = value - value.slough, 
                value.rel = value.diff / value.biom * 100,
                meta = factor(meta)) %>% 
  ggplot() + 
  geom_line(aes(x = timestep, y = value.diff, col = meta)) + 
  facet_wrap(. ~ part, scales = "free_y", ncol = 1) + 
  # geom_hline(yintercept = 0, linetype = 2) + 
  scale_color_viridis_d(name = "", option = "C") +
  labs(x = "Timestep", y = "Cumulative production - slough") + 
  guides(col = guide_legend(nrow = 1)) +
  theme_classic() + theme(legend.position = "bottom")

# required nutrients stable biomass
purrr::map_dfr(result_temp$seafloor, function(i) {
  
  df_sum <- dplyr::select(i, timestep, nutrients_pool, ag_biomass, bg_biomass) %>% 
    dplyr::group_by(timestep) %>%
    dplyr::summarise(nutrients_pool = mean(nutrients_pool),
                     ag_biomass = mean(ag_biomass), bg_biomass = mean(bg_biomass))
  
  list_nutr <- nutr_req(bg_biomass = df_sum$bg_biomass, ag_biomass = df_sum$ag_biomass, 
                        parameters = parameters_list)
  
  data.frame(timestep = df_sum$timestep, 
             bg = df_sum$nutrients_pool - list_nutr$bg, 
             ag = df_sum$nutrients_pool - list_nutr$bg - list_nutr$ag) %>% 
    tidyr::pivot_longer(-timestep)}, .id = "meta") %>% 
  dplyr::mutate(name = factor(name, levels = c("ag", "bg"))) %>% 
  ggplot() + 
  geom_line(aes(x = timestep, y = value, col = meta)) + 
  facet_wrap(. ~ name, scales = "free_y", ncol = 1) +
  # geom_hline(yintercept = 0, linetype = 2) + 
  scale_color_viridis_d(name = "", option = "C") +
  labs(x = "Timestep", y = "Nutrients pool - Stable nutrients") + 
  guides(col = guide_legend(nrow = 1)) +
  theme_classic() + theme(legend.position = "bottom")
