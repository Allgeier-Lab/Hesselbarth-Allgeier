##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

#### Load setup ####

source("05_Various/setup.R")

#### Stable values ####

stable_values <- arrR::get_req_nutrients(bg_biomass = starting_list$bg_biomass,
                                         ag_biomass = starting_list$ag_biomass,
                                         parameters = parameters_list)

starting_list$nutrients_pool <- stable_values$nutrients_pool

starting_list$detritus_pool <- stable_values$detritus_pool

#### Setup HPC ####

df_input <- expand.grid(amplitude_mn = c(0.05, 0.5, 1.0), 
                        enrichment = c(0.5, 1.0, 1.5))

helper_dist <- function(x, y) arrR:::rcpp_closest_reef(x, y, coords_reef = reef_matrix)[[2]]

globals <- c("n", "reef_matrix", "max_i", "starting_list", "parameters_list", "dimensions", 
             "grain", "use_log", "nutrient_input", "freq_mn", "min_per_i", "seagrass_each", 
             "save_each", "helper_dist") 

foo_hpc <- function(amplitude_mn, enrichment) {
  
  library(dplyr)
  
  # setup metaecosystems
  metasyst_temp <- meta.arrR::setup_meta(n = n, max_i = max_i, reef = reef_matrix,
                                         starting_values = starting_list, parameters = parameters_list,
                                         dimensions = dimensions, grain = grain, use_log = use_log, 
                                         verbose = FALSE)
  
  # simulate input
  input_temp <- meta.arrR::simulate_nutr_input(n = n, max_i = max_i, frequency = freq_mn,
                                               input_mn = nutrient_input * enrichment,
                                               amplitude_mn = amplitude_mn, amplitude_sd = 0.0,
                                               phase_mn = 0.0, phase_sd = 0.0, verbose = FALSE)
  
  purrr::map_dfr(list(local = 0, mobile = 8), function(i) {
    
    parameters_list$move_residence_mean <- i
    
    metasyst_temp$fishpop_attr[, 3] <- i
    
    result_temp <- meta.arrR::run_simulation_meta(metasyst = metasyst_temp, parameters = parameters_list,
                                                  nutrients_input = input_temp, movement = "behav",
                                                  max_i = max_i, min_per_i = min_per_i,
                                                  seagrass_each = seagrass_each,
                                                  save_each = save_each, verbose = FALSE)
    
    purrr::map_dfr(result_temp$seafloor, function(j) {
      
      dplyr::filter(j, reef == 0) %>% 
        dplyr::mutate(ttl_biomass = ag_biomass + bg_biomass, ttl_production = ag_production + bg_production, 
                      dist_reef = mapply(helper_dist, x, y), 
                      dist_class = dplyr::case_when(dist_reef <= 3.0 ~ "close",  
                                                    TRUE ~ "far")) %>%
        dplyr::group_by(x, y) %>% 
        dplyr::mutate(ag_production_lag = ag_production - dplyr::lag(ag_production), 
                      bg_production_lag = bg_production - dplyr::lag(bg_production), 
                      ttl_production_lag = ttl_production - dplyr::lag(ttl_production)) %>% 
        dplyr::group_by(timestep, dist_class) %>% 
        dplyr::summarise(ag_production = mean(ag_production), bg_production = mean(bg_production), 
                         ttl_biomass = mean(ttl_biomass), ttl_production = mean(ttl_production),
                         ag_production_lag = mean(ag_production_lag), bg_production_lag = mean(bg_production_lag), 
                         ttl_production_lag = mean(ttl_production_lag), .groups = "drop")}, 
      .id = "meta")}, 
    .id = "movement") %>% 
    dplyr::mutate(amplitude_mn = amplitude_mn, enrichment = enrichment,
                  meta = factor(meta), movement = factor(movement), 
                  dist_class = factor(dist_class, levels = c("close", "far"), 
                                      labels = c("d ≤ 3m", "d > 3 m")))
    
}

#### Submit to HPC #### 

sbatch_pp_dist <- rslurm::slurm_apply(f = foo_hpc, params = df_input, 
                                      global_objects = globals, jobname = "PP_dist",
                                      nodes = nrow(df_input), cpus_per_node = 1, 
                                      slurm_options = list("account" = account, 
                                                           "partition" = "standard",
                                                           "time" = "01:00:00", ## hh:mm::ss
                                                           "mem-per-cpu" = "7G"),
                                      pkgs = c("meta.arrR", "arrR", "dplyr", "purrr"),
                                      rscript_path = rscript_path, sh_template = sh_template, 
                                      submit = FALSE)

#### Collect results ####

suppoRt::rslurm_missing(x = sbatch_pp_dist)

result_pp <- rslurm::get_slurm_out(sbatch_pp_dist, outtype = "table")

suppoRt::save_rds(object = result_pp, filename = "00-result_pp.rds",
                  path = "02_Data/", overwrite = overwrite)

rslurm::cleanup_files(sbatch_pp_dist)

#### Create result figures ####

result_pp <- readr::read_rds("02_Data/00-result_pp.rds")

gg_list <- purrr::map2(df_input$amplitude_mn, df_input$enrichment, function(i, j) {
  
  dplyr::filter(result_pp, amplitude_mn == i, enrichment == j) %>% 
    dplyr::select(-c(amplitude_mn, enrichment)) %>% 
    tidyr::pivot_longer(-c(movement, meta, timestep, dist_class), 
                        names_to = "part") %>% 
    dplyr::group_by(movement, timestep, dist_class, part) %>% 
    dplyr::summarise(mn = mean(value), sd = sd(value), .groups = "drop") %>% 
    dplyr::filter(part %in% c("ag_production_lag", "bg_production_lag", "ttl_production_lag")) %>% 
    ggplot(aes(x = timestep, y = mn, color = movement)) + 
    geom_hline(yintercept = 0, linetype = 2, color = "grey") +
    geom_line() + 
    geom_ribbon(aes(ymin = mn - sd, ymax = mn + sd, fill = movement), col = NA, alpha = 0.25) +
    facet_grid(rows = vars(part), cols = vars(dist_class), scales = "fixed") + 
    scale_color_manual(name = "Cross movement", values = c("#bf6169", "#5e81ac")) +
    scale_fill_manual(name = "Cross movement", values = c("#bf6169", "#5e81ac")) +
    labs(x = "Timestep", y = "Mean PP per sqm", subtitle = paste0("Amplitude mn=", i * 100, "%"), 
         title = paste0("Enrichment=", j * 100, "%")) +
    theme_classic() + theme(panel.border = element_rect(fill = NA, size = 0.25),
                            strip.background = element_rect(fill = NA, size = 0.25))
})
