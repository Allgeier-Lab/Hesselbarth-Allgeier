##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

#### Load setup ####

source("05_Various/setup.R")

source("01_Functions/log_response.R")

#### Adapt parameters ####

# set to null because not needed for arrR
parameters_list$move_residence <- 0.0

#### setup_inputs ####
# create 5 reef cells in center of seafloor
reef_matrix <- matrix(data = c(-1, 0, 0, 1, 1, 0, 0, -1, 0, 0),
                      ncol = 2, byrow = TRUE)

#### setup experiment ####

df_experiment <- expand.grid(enrichment_levels = enrichment_levels, 
                             amplitude_levels = amplitude_levels) %>% 
  dplyr::slice(rep(1:dplyr::n(), each = iterations))

#### Setup function ####

foo <- function(enrichment_levels, amplitude_levels) {
  
  # simulate input
  input_temp <- meta.arrR::sim_nutr_input(n = n, max_i = max_i,
                                          amplitude_mod = amplitude_levels, phase_mod = 0.0,
                                          input_mn = nutrient_input * enrichment_levels,
                                          freq_mn = freq_mn, verbose = FALSE)

  # create list with no fishpop and fishpop
  input_list <- list(0, starting_list$pop_n)
  
  # run model
  result_temp <- purrr::map_dfr(input_list, function(i) {
    
    starting_list$pop_n <- i
    
    metasyst_temp <- meta.arrR::setup_meta(n = n, max_i = max_i, reef = reef_matrix,
                                           starting_values = starting_list, parameters = parameters_list,
                                           dimensions = dimensions, use_log = use_log, 
                                           verbose = FALSE)
    
    # run model
    run_temp <- meta.arrR::run_simulation_meta(metasyst = metasyst_temp, parameters = parameters_list,
                                               nutrients_input = input_temp, movement = "behav",
                                               max_i = max_i, min_per_i = min_per_i,
                                               seagrass_each = seagrass_each, save_each = save_each, 
                                               verbose = FALSE)
    
    result_sum <- meta.arrR::summarize_meta(meta.arrR::filter_meta(x = run_temp, filter = max_i, 
                                                                   verbose = FALSE))
    
    result_sum <- dplyr::left_join(x = result_sum$biomass, y = result_sum$production, 
                                   by = c("meta", "timestep"))
    
    # calculate mean and sum biomass and production
    result_sum <- dplyr::summarise(result_sum, dplyr::across(c(ag_biomass, bg_biomass, ttl_biomass, 
                                                               ag_production, bg_production, ttl_production),
                                                             function(x) c(mean = mean(x), sum = sum(x))))
    
    # add function column
    result_sum <- tibble::add_column(result_sum, fun = c("mean", "sum"), .before = "ag_biomass") 
    
    # reshape to longer format
    seafloor_temp <- tidyr::pivot_longer(result_sum, -fun) 
    
    # add treatment levels
    dplyr::mutate(seafloor_temp, enrichment = enrichment_levels, amplitude = amplitude_levels, 
                  fishpop = ifelse(test = unique(run_temp$starting_values$pop_n) == 0, 
                                   yes = "nofish", no = "local"))
    })
  
  # reshape wider using nofish and local
  result_temp <- tidyr::pivot_wider(result_temp, names_from = fishpop, values_from = value)
  
  return(result_temp)
  
}

#### Submit to HPC ####

globals <- c("n", "max_i", "nutrient_input", "freq_mn", 
             "starting_list", 
             "reef_matrix", "parameters_list", "dimensions", "use_log", 
             "min_per_i", "seagrass_each", "save_each")

# create .sh script
sbatch_pp <- rslurm::slurm_apply(f = foo, params = df_experiment,
                                 global_objects = globals, jobname = "PP_nofish_local",
                                 nodes = nrow(df_experiment), cpus_per_node = 1, 
                                 slurm_options = list("account" = account, 
                                                      "partition" = "standard",
                                                      "time" = "01:00:00", ## hh:mm::ss
                                                      "mem-per-cpu" = "7G"),
                                 pkgs = c("arrR", "purrr", "dplyr", "tidyr", "tibble"),
                                 rscript_path = rscript_path, sh_template = sh_template, 
                                 submit = FALSE)

# check results
suppoRt::rslurm_missing(sbatch_pp)

# get results as list
pp_nofish_local <- rslurm::get_slurm_out(sbatch_pp, outtype = "table")

# save results to disk
suppoRt::save_rds(object = pp_nofish_local, file = "02_Data/00_pp_nofish_local.rds", 
                  overwrite = overwrite)

# delete .sh scripts
rslurm::cleanup_files(sbatch_pp)

#### Settings ggplot ####

parts_list <- list(ag = c("ag_biomass", "ag_production"), bg = c("bg_biomass", "bg_production"), 
                   ttl = c("ttl_biomass", "ttl_production"))

legend_postion <- c("none", "none", "bottom")

x_labels <- list("", "", "Enrichment treatment")

x_axis <- list(element_blank(), element_blank(), NULL)

col_palette <- c("#5ABCD6", "#FAD510", "#F22301")

#### Results treatments (Rel difference)  ####

pp_nofish_local <- readr::read_rds("02_Data/00_pp_nofish_local.rds") %>%
  dplyr::mutate(reldiff = ((local - nofish) / nofish) * 100) %>%
  dplyr::mutate(fun = factor(fun), name = factor(name),
                enrichment = factor(enrichment, ordered = TRUE),
                amplitude = factor(amplitude, ordered = TRUE))

y_labels <- list("", expression(paste(Delta, "no fish, local fish [%]")), "")

gg_list <- map(seq_along(parts_list), function(i) {
  
  dplyr::filter(pp_nofish_local, fun == "sum", name %in% parts_list[[i]]) %>% 
    ggplot() +
    geom_hline(yintercept = 0, linetype = 2, color = "grey") +
    geom_boxplot(aes(x = enrichment, y = reldiff, fill = amplitude),
                 position = position_dodge(), alpha = 0.75) +
    geom_point(aes(x = enrichment, y = reldiff, col = amplitude),
               position = position_jitterdodge(), alpha = 0.25) +
    facet_wrap(. ~ name, nrow = 1) +
    labs(x = x_labels[[i]], y = y_labels[[i]]) +
    scale_fill_manual(name = "Amplitude treatment", values = col_palette) +
    scale_color_manual(name = "Amplitude treatment", values = col_palette) +
    theme_classic() + theme(legend.position = legend_postion[[i]], 
                            axis.text.x = x_axis[[i]])
  
})

cowplot::plot_grid(plotlist = gg_list, nrow = 3, rel_heights = c(0.75, 0.75, 1))

# dplyr::group_by(pp_nofish_local, fun, name, enrichment, amplitude) %>%
#   dplyr::summarise(reldiff_mn = mean(reldiff), reldiff_sd = sd(reldiff),
#                    .groups = "drop") %>%
#   dplyr::arrange(desc(fun), name, enrichment, amplitude) %>%
#   dplyr::filter(fun == "sum")

#### Results treatments (Response ratios)  ####

pp_nofish_local <- readr::read_rds("02_Data/00_pp_nofish_local.rds") %>% 
  dplyr::group_by(fun, name, enrichment, amplitude) %>% 
  dplyr::group_split() %>% 
  purrr::map_dfr(function(i) {
    
    bootstrap <- boot::boot(data = tibble::tibble(ctrl = i$nofish, 
                                                  trtm = i$local),
                            statistic = log_response, R = 10000)
    
    bootstrap_ci <- boot::boot.ci(bootstrap, type = "norm", conf = 0.95)
    
    tibble(fun = unique(i$fun), name = unique(i$name), 
           enrichment = unique(i$enrichment), amplitude = unique(i$amplitude),
           mean = mean(bootstrap$t[, 1]), lo = bootstrap_ci$normal[2], hi = bootstrap_ci$normal[3])}) %>% 
  dplyr::mutate(fun = factor(fun), name = factor(name),
                enrichment = factor(enrichment, ordered = TRUE),
                amplitude = factor(amplitude, ordered = TRUE))

y_labels <- list("", "log(RR) no fish, local fish", "")

gg_list <- map(seq_along(parts_list), function(i) {
  
  dplyr::filter(pp_nofish_local, fun == "sum", name %in% parts_list[[i]]) %>% 
    ggplot() +
    geom_hline(yintercept = 0, linetype = 2, color = "grey") +
    geom_point(aes(x = enrichment, y = mean, col = amplitude), 
               position = position_dodge(width = 0.5), alpha = 0.25) + 
    geom_errorbar(aes(x = enrichment, ymin = lo, ymax = hi, col = amplitude),  
                  position = position_dodge(width = 0.5), width = 0.25) +
    facet_wrap(. ~ name, ncol = 2) +
    labs(x = x_labels[[i]], y = y_labels[[i]]) +
    scale_fill_manual(name = "Amplitude treatment", values = col_palette) +
    scale_color_manual(name = "Amplitude treatment", values = col_palette) +
    theme_classic() + theme(legend.position = legend_postion[[i]], 
                            axis.text.x = x_axis[[i]])
  
})

cowplot::plot_grid(plotlist = gg_list, nrow = 3, rel_heights = c(0.75, 0.75, 1))
