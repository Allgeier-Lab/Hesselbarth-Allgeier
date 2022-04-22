##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

#### Load setup ####

source("05_Various/setup.R")

#### Adapt parameters ####

parameters_list$move_residence_mean <- 0.0

#### Stable values #### 

stable_values <- arrR::get_req_nutrients(bg_biomass = starting_list$bg_biomass,
                                         ag_biomass = starting_list$ag_biomass,
                                         parameters = parameters_list)

starting_list$nutrients_pool <- stable_values$nutrients_pool

starting_list$detritus_pool <- stable_values$detritus_pool

#### Setup experiment ####

# create random sd for each iteration of amplitude levels and stochastic treatments
rndm_sd <- runif(n = iterations * length(amplitude_levels) * 3, min = 0.0, max = 1.0)

# create data.frame with all combinations
experiment_df <- tibble::tibble(ampl_mn = rep(x = amplitude_levels, each = iterations * 3), 
                                rndm_sd = rndm_sd, 
                                stochastic = rep(x = rep(x = c("amplitude", "phase", "both"), 
                                                         each = iterations), times = 3))

# setup HPC function
foo_hpc <- function(ampl_mn, rndm_sd, stochastic) {
  
  # setup metaecosystems
  metasyst_temp <- meta.arrR::setup_meta(n = n, max_i = max_i, reef = reef_matrix,
                                         starting_values = starting_list, parameters = parameters_list,
                                         dimensions = dimensions, grain = grain, use_log = use_log, 
                                         verbose = FALSE)
  
  # get sd values depending on what to vary
  amplitude_sd <- ifelse(test = stochastic == "phase", yes = 0, no = rndm_sd)
    
  phase_sd <- ifelse(test = stochastic == "amplitude", yes = 0, no = rndm_sd)

  # simulate nutrient input
  input_temp <-  meta.arrR::simulate_nutr_input(n = n, max_i = max_i, frequency = freq_mn, 
                                                input_mn = nutrient_input, 
                                                amplitude_mn = ampl_mn, amplitude_sd = amplitude_sd, 
                                                phase_mn = 0.0, phase_sd = phase_sd, 
                                                verbose = FALSE)

  # run model
  result_temp <- meta.arrR::run_simulation_meta(metasyst = metasyst_temp, parameters = parameters_list,
                                                nutrients_input = input_temp, movement = "behav",
                                                max_i = max_i, min_per_i = min_per_i,
                                                seagrass_each = seagrass_each,
                                                save_each = save_each, verbose = FALSE)
  
  # calc cv
  cv <- meta.arrR::calc_variability(x = result_temp, lag = c(FALSE, TRUE))
  
  # combine to result data.frame
  dplyr::mutate(dplyr::bind_rows(cv), ampl_mn = ampl_mn, rndm_sd = rndm_sd, stochastic = stochastic)
                
}

#### Submit HPC

globals <- c("n", "max_i", "reef_matrix", "starting_list", "parameters_list", "dimensions", 
             "grain", "use_log", "freq_mn", "nutrient_input", "min_per_i", "seagrass_each", 
             "save_each") 

sbatch_cv <- rslurm::slurm_apply(f = foo_hpc, params = experiment_df, 
                                 global_objects = globals, jobname = "pp_abiotic_cv",
                                 nodes = nrow(experiment_df), cpus_per_node = 1, 
                                 slurm_options = list("account" = account, 
                                                      "partition" = "standard",
                                                      "time" = "01:00:00", ## hh:mm::ss
                                                      "mem-per-cpu" = "7G"),
                                 pkgs = c("dplyr", "meta.arrR"),
                                 rscript_path = rscript_path, sh_template = sh_template, 
                                 submit = FALSE)

#### Collect results #### 

suppoRt::rslurm_missing(x = sbatch_cv)

cv_result <- rslurm::get_slurm_out(sbatch_cv, outtype = "table")

suppoRt::save_rds(object = cv_result, filename = "01-abiotic-variability.rds",
                  path = "02_Data/", overwrite = overwrite)

rslurm::cleanup_files(sbatch_cv)

#### Results ####

cv_result <- readr::read_rds("02_Data/01-abiotic-variability.rds") %>% 
  dplyr::mutate(part = factor(part), measure = factor(measure), 
                ampl_mn = factor(ampl_mn, ordered = TRUE, 
                                 labels = c("Low amplitude mean", "Medium amplitude mean", "High amplitude mean")),
                stochastic = factor(stochastic, levels = c("amplitude", "phase", "both"), 
                                    labels = c("Amplitude variability", "Phase variability", "Simultaneous variability"))) %>% 
  tibble::tibble()

dplyr::filter(cv_result, measure %in% c("alpha", "gamma"), 
              part %in% c("ag_production", "bg_production", "ttl_production"),
              stochastic != "Simultaneous variability") %>% 
  ggplot(aes(x = rndm_sd, y = value, linetype = measure, color = ampl_mn)) + 
  geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
  geom_point(alpha = 0.25, shape = 19, size = 1.5) + 
  geom_smooth(method = "loess", se = FALSE, size = 0.75) +
  facet_grid(rows = vars(part), cols = vars(stochastic), scales = "free_y") + 
  scale_x_continuous(breaks = seq(from = 0, to = 1.0, length.out = 5), limits = c(0, 1)) +
  scale_color_manual(name = "", values = c("#007c2aff", "#ffcc00ff", "#2f62a7ff")) +
  scale_linetype_manual(name = "", values = c(1, 2), 
                        labels = c(expression(paste(alpha, " scale")), 
                                   expression(paste(gamma, " scale")))) +
  labs(x = "Abiotic variability", y = expression("CV"["primary production"])) + 
  guides(linetype = guide_legend(order = 1, keywidth = unit(10, "mm"), 
                                 override.aes = list(color = "black"), nrow = 2), 
         color = guide_legend(order = 2, nrow = 3)) +
  theme_classic() + theme(legend.position = "bottom", text = element_text(family = "Georgia"))

dplyr::filter(cv_result, measure %in% c("alpha", "gamma"), 
              part %in% c("ag_production", "bg_production", "ttl_production"), 
              stochastic != "Simultaneous variability") %>% 
  tidyr::pivot_wider(names_from = measure, values_from = value) %>% 
  dplyr::mutate(beta = alpha / gamma) %>% 
  ggplot(aes(x = alpha, y = gamma, color = ampl_mn)) + 
  # geom_polygon(data = data.frame(x = c(0.0, 0.5, 0.5), y = c(0.0, 0.0, 0.5)), 
  #              aes(x = x, y = y), fill = "#a8bd90", col = NA, alpha = 0.65) +
  # geom_polygon(data = data.frame(x = c(0.0, 0.0, 0.5), y = c(0.0, 0.5, 0.5)), 
  #              aes(x = x, y = y), fill = "#b3666b", col = NA, alpha = 0.65) +
  annotate(geom = "text", x = 0.25, y = 0.275, angle = 45, 
           color = "black", label = "paste(beta, '=1')", parse = TRUE) +
  annotate(geom = "text", x = 0.1, y = 0.45, color = "#b3666b", size = 2.5,
           label = "negative PE") +
  annotate(geom = "text", x = 0.4, y = 0.05, color = "#a8bd90", size = 2.5,
           label = "positive PE") +
  geom_point(shape = 19, alpha = 0.5) +
  geom_abline(slope = 1, linetype = 2, color = "grey") +
  facet_grid(rows = vars(part), cols = vars(stochastic)) + 
  scale_x_continuous(limits = c(0, 0.5)) + scale_y_continuous(limits = c(0, 0.5)) +
  scale_color_manual(name = "", values = c("#007c2a", "#ffcc00", "#2f62a7")) + 
  labs(x = expression(paste(alpha, " CV"["primary production"])),
       y = expression(paste(gamma, " CV"["primary production"]))) +
  guides(size = "none", color = guide_legend(nrow = 1)) +
  coord_equal() +
  theme_classic() + theme(legend.position = "bottom")
