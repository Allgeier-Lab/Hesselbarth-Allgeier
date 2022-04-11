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

parameters_list$move_residence <- 0.0

parameters_list$move_residence_var <- 0.0

#### Stable values ####

stable_values <- arrR::get_stable_values(bg_biomass = starting_list$bg_biomass,
                                         ag_biomass = starting_list$ag_biomass,
                                         parameters = parameters_list)

starting_list$nutrients_pool <- stable_values$nutrients_pool

starting_list$detritus_pool <- stable_values$detritus_pool

#### Simulate input ####

# setup enrichment levels
enrichment_levels <- rep(c(low = 0.75, medium = 1.0, high = 1.25), each = iterations)

variability <- runif(n = length(enrichment_levels), min = 0.0, max = 1.0)

df_experiment <- data.frame(variability = variability, enrichment = enrichment_levels)

#### Setup HPC function ####

globals <- list(n = n, reef_matrix = reef_matrix, max_i = max_i, starting_list = starting_list, 
                parameters_list = parameters_list, dimensions = dimensions, 
                grain = grain, use_log = use_log, input_mn = stable_values$nutrients_input, 
                freq_mn = freq_mn, min_per_i = min_per_i, seagrass_each = seagrass_each, 
                save_each = save_each) 

foo <- function(variability, enrichment) {
  
  # setup metaecosystems
  metasyst_temp <- meta.arrR::setup_meta(n = globals$n, max_i = globals$max_i, reef = globals$reef_matrix,
                                         starting_values = globals$starting_list,
                                         parameters = globals$parameters_list,
                                         dimensions = globals$dimensions, grain = globals$grain,
                                         use_log = globals$use_log, verbose = FALSE)
  
  # simulate input
  input_temp <- meta.arrR::sim_nutr_input(n = globals$n, max_i = globals$max_i,
                                          variability = variability, input_mn = globals$input_mn * enrichment, 
                                          freq_mn = globals$freq_mn)
  
  # run model
  result_temp <- meta.arrR::run_simulation_meta(metasyst = metasyst_temp, parameters = globals$parameters_list,
                                                nutrients_input = input_temp, movement = "behav",
                                                max_i = globals$max_i, min_per_i = globals$min_per_i,
                                                seagrass_each = globals$seagrass_each,
                                                save_each = globals$save_each, verbose = FALSE)
  
  # filter only second half of timesteps
  result_temp <- meta.arrR::filter_meta(x = result_temp, filter = c(globals$max_i / 2,
                                                                    globals$max_i), 
                                        reset = TRUE, verbose = FALSE)
  
  # calc CV
  cv_temp <- dplyr::bind_rows(meta.arrR::calc_variability(x = result_temp, lag = c(FALSE, TRUE), 
                                                          verbose = FALSE))
  
  # add type col
  cv_temp <- dplyr::bind_cols(type = "cv", cv_temp)
  
  # get the total sum of each timestep and metaecosystem
  result_sum <- dplyr::bind_rows(meta.arrR::summarize_meta(result = result_temp))
  
  # filter only final timestep
  result_sum <- dplyr::filter(result_sum, timestep == max(timestep))
  
  # calculate mean (alpha) and total (gamma) biomass values
  result_sum <- apply(X = result_sum[, -c(1, 2)] , MARGIN = 2, 
                      FUN = function(x) c(alpha = mean(x, na.rm = TRUE), 
                                          gamma = sum(x, na.rm = TRUE)))
  
  # reshape to long format
  result_sum <- tidyr::pivot_longer(data.frame(type = "absolute", measure = row.names(result_sum), 
                                               result_sum), -c(type, measure), 
                                    names_to = "part")
  
  # rearrange cols
  result_sum <- dplyr::select(result_sum, type, part, measure, value)
  
  # combine to one df
  dplyr::bind_cols(enrichment = enrichment, variability = variability, 
                   dplyr::arrange(dplyr::bind_rows(cv_temp, result_sum), 
                                  type, part, measure))
  
}

#### Submit to HPC #### 

sbatch_pe_prod <- rslurm::slurm_apply(f = foo, params = df_experiment, 
                                      global_objects = "globals", jobname = "PE_Prod_local",
                                      nodes = nrow(df_experiment), cpus_per_node = 1, 
                                      slurm_options = list("account" = account, 
                                                           "partition" = "standard",
                                                           "time" = "02:00:00", ## hh:mm::ss
                                                           "mem-per-cpu" = "7G", 
                                                           "exclude" = exclude_nodes),
                                      pkgs = c("dplyr", "meta.arrR", "tidyr"),
                                      rscript_path = rscript_path, sh_template = sh_template, 
                                      submit = FALSE)

#### Collect results ####

suppoRt::rslurm_missing(sbatch_pe_prod)

df_pe_prod <- rslurm::get_slurm_out(sbatch_pe_prod, outtype = "table")

suppoRt::save_rds(object = df_pe_prod, filename = "02_PE-Prod_Enrich-local.rds", 
                  path = "02_Data/", overwrite = overwrite)

rslurm::cleanup_files(sbatch_pe_prod)

#### Load data ####

df_pe_prod <- readRDS("02_Data/02_PE-Prod_Enrich-local.rds") %>% 
  dplyr::filter(c(type == "cv" & measure == "beta") | c(type == "absolute" & measure == "gamma")) %>% 
  dplyr::select(-type) %>%
  tidyr::pivot_wider(names_from = "measure", values_from = value, values_fn = list) %>% 
  tidyr::unnest(cols = tidyselect::everything()) %>% 
  dplyr::mutate(enrichment = factor(enrichment, levels = c(0.75, 1.0, 1.25), labels = c("low", "medium", "high")),
                variability = as.numeric(variability),
                part = factor(part, levels = c("ag_biomass", "ag_production", "ttl_biomass",
                                               "bg_biomass", "bg_production", "ttl_production")))

#### Create CV ggplot ####

labels_facet <- list(c(low = "Low enrichment", medium = "Medium enrichment", high = "High enrichment"), 
                     c(low = "", medium = "", high = ""), 
                     c(low = "", medium = "", high = ""))

x_axis <- c(" ", " ", "Portfolio effect (beta)")

y_axis <- c(" ", "log Abs value per sqm", " ")

gg_results <- purrr::map(c("production", "biomass"), function(i) { 
  
  # create parts to loop through
  parts <- paste0(c("ag_", "bg_", "ttl_"), i)
  
  # create names for plot labelling
  names(parts) <- paste(c("Aboveground", "Belowground", "Total"), i)
  
  gg_temp <- purrr::map(seq_along(parts), function(i){
    
    dplyr::filter(df_pe_prod, part %in% parts[[i]]) %>% 
      dplyr::mutate(gamma = log(gamma / prod(dimensions))) %>% 
      ggplot(aes(x = beta, y = gamma)) +
      geom_point(shape = 1) + 
      geom_smooth(method = lm, se = TRUE, linetype = 2, col = "black", formula = 'y ~ x') +
      stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`;`~")),
               p.accuracy = 0.01, r.accuracy = 0.001) +
      facet_wrap(. ~ enrichment, ncol = 3, nrow = 1, scales = "fixed",
                 labeller = labeller(enrichment = labels_facet[[i]])) + 
      labs(y = y_axis[[i]], x = x_axis[[i]], subtitle = names(parts)[i]) +
      theme_classic(base_size = base_size) + 
      theme(legend.position = "bottom", strip.background = element_blank(), 
            strip.text = element_text(hjust = 0))
    
  })
  
  cowplot::plot_grid(plotlist = gg_temp, ncol = 1, nrow = 3)

})

#### Save ggplot ####

# set name of output level
names(gg_results) <- c("production", "biomass")

#### Save ggplot ####

# # loop through output level
# purrr::walk(seq_along(gg_results), function(i) {
#   
#   # create file name
#   filename_temp <- paste0("02_pe_prod_", names(gg_results)[[i]], "_local.png")
#     
#   # save ggplot
#   suppoRt::save_ggplot(plot = gg_results[[i]], filename = filename_temp,
#                        path = "04_Figures", width = height, height = width, dpi = dpi, 
#                        units = units, overwrite = overwrite)
# })
