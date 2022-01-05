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

default_starting$pop_n <- 0

default_parameters$nutrients_diffusion <- 0.0
default_parameters$detritus_diffusion <- 0.0
default_parameters$detritus_fish_diffusion <- 0.0

# default_parameters$seagrass_thres <- 1/3

#### Stable values ####

stable_values <- arrR::get_stable_values(starting_values = default_starting,
                                         parameters = default_parameters)

default_starting$nutrients_pool <- stable_values$nutrients_pool

default_starting$detritus_pool <- stable_values$detritus_pool

#### Simulate input ####

# number of iterations
itr <- 50

# setup enrichment levels
enrichment_levels <- rep(c(low = 0.5, medium = 0.75, high = 1.0), each = itr)

variability <- runif(n = length(enrichment_levels), min = 0.0, max = 1.0)

df_experiment <- data.frame(variability = variability, enrichment = enrichment_levels)

#### Setup HPC function ####

globals <- list(n = n, max_i = max_i, default_starting = default_starting, 
                default_parameters = default_parameters, dimensions = dimensions, 
                grain = grain, input_mn = stable_values$nutr_input, freq_mn = freq_mn,
                min_per_i = min_per_i, seagrass_each = seagrass_each, save_each = save_each) 

foo <- function(variability, enrichment) {
  
  # setup metaecosystems
  metasyst_temp <- meta.arrR::setup_meta(n = globals$n, max_i = globals$max_i,
                                         starting_values = globals$default_starting,
                                         parameters = globals$default_parameters,
                                         dimensions = globals$dimensions, grain = globals$grain,
                                         reefs = NULL, verbose = FALSE)
  
  # simulate input
  input_temp <- meta.arrR::sim_nutr_input(n = globals$n, max_i = globals$max_i,
                                          variability = variability,
                                          input_mn = globals$input_mn * enrichment, 
                                          freq_mn = globals$freq_mn)
  
  # run model
  result_temp <- meta.arrR::run_meta(metasyst = metasyst_temp, nutr_input = input_temp,
                                     parameters = globals$default_parameters,
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
                                      global_objects = "globals", jobname = "PE_Prod",
                                      nodes = nrow(df_experiment), cpus_per_node = 1, 
                                      slurm_options = list("account" = "jeallg1", 
                                                           "partition" = "standard",
                                                           "time" = "02:00:00", ## hh:mm::ss
                                                           "mem-per-cpu" = "7G", 
                                                           "exclude" = exclude_nodes),
                                      pkgs = c("dplyr", "meta.arrR", "tidyr"),
                                      rscript_path = rscript_path, sh_template = sh_template, 
                                      submit = FALSE)

#### Collect results ####

df_pe_prod <- rslurm::get_slurm_out(sbatch_pe_prod, outtype = "table")

suppoRt::save_rds(object = df_pe_prod, filename = "02_PE-Prod_Enrich.rds", 
                  path = "02_Data/", overwrite = overwrite)

rslurm::cleanup_files(sbatch_pe_prod)

#### Load data ####

df_pe_prod <- readRDS("02_Data/02_PE-Prod_Enrich.rds") %>% 
  dplyr::filter(c(type == "cv" & measure == "beta") | c(type == "absolute" & measure == "gamma")) %>% 
  dplyr::select(-type) %>%
  tidyr::pivot_wider(names_from = "measure", values_from = value, values_fn = list) %>% 
  tidyr::unnest(cols = tidyselect::everything()) %>% 
  dplyr::mutate(enrichment = factor(enrichment, levels = c(0.5, 0.75, 1.0), labels = c("low", "medium", "high")),
                variability = as.numeric(variability),
                part = factor(part, levels = c("ag_biomass", "ag_production", "ttl_biomass",
                                               "bg_biomass", "bg_production", "ttl_production")))

#### Create CV ggplot ####

# col_palette_part <- c("#586F7E", "#168B98", "#ED5B66")
# 
# col_palette_enrich <- c("#C7247B", "#A5E100", "#172969")
# 
# col_palette <- c("#5ABCD6", "#FAD510", "#F22301")

# create some ggplot settings/objects
parts <- c(Aboveground = "ag_production", Belowground = "bg_production", Total = "ttl_production")

labels_facet <- list(c(low = "Low enrichment", medium = "Medium enrichment", high = "High enrichment"), 
                     c(low = "", medium = "", high = ""), 
                     c(low = "", medium = "", high = ""))

x_axis <- c(" ", " ", "Portfolio effect (beta)")

y_axis <- c(" ", "Absolute production per sqm", " ")

gg_pe_prod <- purrr::map(seq_along(parts), function(i){
  
  dplyr::filter(df_pe_prod, part == parts[i]) %>% 
    dplyr::mutate(gamma = gamma / 10000) %>% 
    ggplot(aes(x = beta, y = gamma)) +
    geom_point(pch = 1) + 
    geom_smooth(method = lm, se = TRUE, col = "black") + 
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`;`~")),
             p.accuracy = 0.01, r.accuracy = 0.001) +
    facet_wrap(. ~ enrichment, ncol = 3, nrow = 1, scales = "free_x",
               labeller = labeller(enrichment = labels_facet[[i]])) + 
    # expand_limits(x = 1) + 
    scale_color_manual(name = "", values = col_palette_part) +
    labs(y = y_axis[[i]], x = x_axis[[i]], subtitle = names(parts)[i]) +
    theme_classic(base_size = base_size) + 
    theme(legend.position = "none", strip.background = element_blank(), 
          strip.text = element_text(hjust = 0))
  
})

gg_pe_prod <- cowplot::plot_grid(plotlist = gg_pe_prod, ncol = 1, nrow = 3)

dplyr::mutate(df_pe_prod, enrichment = factor(enrichment, 
                                              labels = c("Low enrichment", 
                                                         "Medium enrichment", 
                                                         "High enrichment"))) %>% 

gg_var_density <- dplyr::mutate(df_pe_prod, enrichment = factor(enrichment,
                                                                labels = c("Low enrichment", 
                                                                           "Medium enrichment", 
                                                                           "High enrichment"))) %>% 
  ggplot() +
  geom_density(aes(x = variability, fill = enrichment, col = enrichment), alpha = 0.25) + 
  labs(x = "Input variability", y = "Density") +
  scale_color_manual(name = "", values = c("#C7247B", "#A5E100", "#172969")) +
  scale_fill_manual(name = "", values = c("#C7247B", "#A5E100", "#172969")) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

#### Save ggplot ####

suppoRt::save_ggplot(plot = gg_pe_prod, filename = "02_pe_prod.png",
                     path = "04_Figures", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)

suppoRt::save_ggplot(plot = gg_var_density, filename = "02_var_density.png",
                     path = "04_Figures", width = height, height = width, dpi = dpi, 
                     units = units, overwrite = overwrite)
