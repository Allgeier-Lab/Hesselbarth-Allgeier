##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

#### Load setup ####

source("01_Functions/setup.R")
# library(future)
# library(furrr)

source("01_Functions/import-cv.R")
source("01_Functions/import-mortality.R")

mem_per_cpu <- "7G"
time <- "05:00:00" # hh:mm:ss

#### Setup experiment ####

# experiment_df <- readRDS("02_Data/experiment-parameters.rds")

pop_n <- c(8, 16, 32, 64, 128)

nutrient_input_fish <- readRDS("02_Data/nutrient-input-fish.rds") |> 
  dplyr::group_by(pop_n) |> 
  dplyr::summarise(excretion_cell = mean(excretion_cell)) |> 
  dplyr::mutate(pop_n = factor(pop_n, ordered = TRUE))

nutrient_levels <- c(low = nutrient_input_fish$excretion_cell[5] * 1.5, medium = 0.000121005, high = 0.000191781)

experiment_df <- expand.grid(pop_n = pop_n, nutrient_input = nutrient_levels) |> 
  dplyr::slice(rep(1:dplyr::n(), each = iterations))

#### Init HPC function ####

foo_hpc <- function(pop_n, nutrient_input) {

  starting_values_list$pop_n <- pop_n

  # # update move meta_sd parameters
  # parameters_list$move_meta_sd <- biotic

  # setup metaecosystems
  metasyst_temp <- meta.arrR::setup_meta(n = n, max_i = max_i, reef = reef_matrix,
                                         starting_values = starting_values_list, parameters = parameters_list,
                                         dimensions = dimensions, grain = grain,
                                         use_log = FALSE, verbose = FALSE)

  # simulate nutrient input
  input_temp <- meta.arrR::simulate_nutrient_noise(n = n, max_i = max_i, input_mn = nutrient_input,
                                                   verbose = FALSE)

  # run model
  result_temp <- meta.arrR::run_simulation_meta(metasyst = metasyst_temp, parameters = parameters_list,
                                                nutrients_input = input_temp, movement = "behav",
                                                torus_diffusion = TRUE, max_i = max_i, min_per_i = min_per_i,
                                                seagrass_each = seagrass_each, save_each = save_each,
                                                verbose = FALSE)  |>
    meta.arrR::filter_meta(filter = c((max_i / years) * years_filter, max_i),
                           reset = TRUE, verbose = FALSE)

  # calc mortality
  mortality <- dplyr::bind_rows(result_temp$fishpop) |>
    dplyr::filter(timestep == max_i) |>
    dplyr::mutate(died_total = died_consumption + died_background) |>
    dplyr::select(id, died_consumption, died_background, died_total) |>
    dplyr::arrange(id)

  prod <- meta.arrR::summarize_meta(result = result_temp, biomass = FALSE, production = TRUE,
                                    fun = function(x, ...) {mean(x, ...) / ((save_each * 120) / 60 / 24)},
                                    lag = c(NA, TRUE))[["production"]] |>
    dplyr::filter(timestep != min(timestep)) |>
    tidyr::pivot_longer(-c(meta, timestep), names_to = "part") |>
    dplyr::group_by(timestep, part) |>
    dplyr::summarise(alpha = mean(value), gamma = sum(value), .groups = "drop") |>
    tidyr::pivot_longer(-c(timestep, part), names_to = "measure") |>
    dplyr::group_by(part, measure) |>
    dplyr::summarise(value = mean(value), .groups = "drop")

  # calc cv
  cv <- meta.arrR::calc_variability(x = result_temp, biomass = FALSE, production = TRUE,
                                    fun = function(x, ...) {mean(x, ...) / ((save_each * 120) / 60 / 24)},
                                    lag = c(NA, TRUE))[["production"]]

  # combine to result data.frame and list
  list(fishpop_init = dplyr::mutate(dplyr::bind_rows(metasyst_temp$fishpop), pop_n = pop_n, nutrient_input = nutrient_input),
       mortality = dplyr::mutate(mortality, pop_n = pop_n, nutrient_input = nutrient_input),
       prod = dplyr::mutate(prod, pop_n = pop_n, nutrient_input = nutrient_input),
       cv = dplyr::mutate(dplyr::bind_rows(cv), pop_n = pop_n, nutrient_input = nutrient_input)
  )
}

# #### Run locally using future ####
# 
# future::plan(multisession, workers = 12)
# 
# results_null <- future_map(1:nrow(experiment_df), function(i){
#   foo_hpc(pop_n = experiment_df$pop_n[i], nutrient_input = experiment_df$nutrient_input[i])
# }, .options = furrr_options(seed = TRUE))
# 
# suppoRt::save_rds(object = results_null, path = "02_Data/", filename = "result-null.rds", 
#                   overwrite = FALSE)

#### Submit HPC ####

globals <- c("n", "max_i", "reef_matrix", "starting_values_list", "parameters_list", "dimensions", "grain", # setup_meta
             "min_per_i", "seagrass_each", "save_each", # run_simulation_meta
             "years", "years_filter") # filter_meta

sbatch_noise <- rslurm::slurm_apply(f = foo_hpc, params = experiment_df,
                                    global_objects = globals, jobname = "null_model",
                                    nodes = nrow(experiment_df), cpus_per_node = 1,
                                    slurm_options = list("account" = account,
                                                         "partition" = "standard",
                                                         "time" = time,
                                                         "mem-per-cpu" = mem_per_cpu),
                                    pkgs = c("arrR", "dplyr", "meta.arrR", "purrr", "tidyr"),
                                    rscript_path = rscript_path, submit = FALSE)

#### Collect results ####

suppoRt::rslurm_missing(x = sbatch_noise)

cv_noise <- rslurm::get_slurm_out(sbatch_noise, outtype = "raw")

suppoRt::save_rds(object = cv_noise, path = "02_Data/", filename = "result-null.rds",
                  overwrite = FALSE)

rslurm::cleanup_files(sbatch_noise)

#### Load data ### 

results_null <- import_cv_null("02_Data/result-null.rds")

mortality_df <- import_mortality_null(path = "02_Data/result-null.rds") |> 
  dplyr::group_by(row_id, pop_n, nutrient_input) |>
  dplyr::summarise(died_total = mean(died_total), .groups = "drop")

#### Filter and process data ####

results_final_df <- dplyr::left_join(x = results_null, y = mortality_df, 
                                     by = c("row_id", "pop_n", "nutrient_input")) |> 
  dplyr::mutate(include = dplyr::case_when(died_total > threshold_mort ~ "no", TRUE ~ "yes")) |> 
  dplyr::filter(include == "yes") |> 
  dplyr::mutate(part = factor(part, levels = c("ag_production", "bg_production", "ttl_production"),
                              labels = c("ag_production" = "Aboveground", "bg_production" = "Belowground", 
                                         "ttl_production" = "Total")), 
                measure = factor(measure, levels = c("alpha", "gamma", "beta", "synchrony"), 
                                 labels = c("Local", "Meta-Ecosys.", "Portfolio effect", "Synchrony"))) |> 
  dplyr::group_by(part, measure, nutrient_input) |> #pop_n
  dplyr::summarise(mn = mean(value.cv), sd = sd(value.cv)) |> 
  dplyr::mutate(facet = dplyr::case_when(measure %in% c("Local", "Meta-Ecosys.") ~ "left", 
                                         TRUE ~ "right"))

#### Create ggplot #### 

base_size <- 13.0
color_enrich <- c("low" = "#009E73", "medium" = "#F0E442", "high" = "#CC79A7")
w <- 0.5

gg_null_list <- purrr::map(c("Aboveground", "Total"), function(i) { 
  
  lgd_post <-  c(0.9, 0.9)
  label_greek <- data.frame(facet = c("left", "right"), label = c("i)", "ii)"),
                            x = c(1, 1), y = c(0.15, 3))

  if(i == "Total") {
    lgd_post <- "none"
    label_greek <- data.frame(facet = c("left", "right"), label = c("iii)", "iv)"),
                              x = c(1, 1), y = c(0.0115, 2.55))
  }
  
  dplyr::filter(results_final_df, part == i, measure %in% c("Local", "Meta-Ecosys.", "Portfolio effect")) |>
    ggplot(aes(x = measure)) +
    
    # adding geoms
    # mean stability values bars
    geom_col(aes(y = mn, fill = nutrient_input), position = position_dodge(w * 0.75),
             width = w * 0.70, color = NA) +
    
    geom_errorbar(aes(ymax = mn + sd, ymin = mn - sd, color = nutrient_input),
                  position = position_dodge(w * 0.75), width = 0.0, linewidth = 0.5) +
    
    # adding facet labels
    geom_text(data = label_greek, aes(x = x , y = y, label = label), size = 5.5, 
              nudge_x = -0.35) +
    
    # scales and labs
    scale_color_manual(name = "", values = color_enrich) +
    scale_fill_manual(name = "", values = color_enrich) +
    labs(x = "", y = "") +
    
    # adding box
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
    annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
    annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
    
    # facet_wrap by part
    facet_wrap(. ~ facet, scales = "free") +
    
    # themeing
    theme_classic(base_size = base_size) +
    theme(legend.position = lgd_post, legend.background = element_blank(),
          legend.text = element_text(size = 6.5),
          axis.line = element_blank(), strip.background = element_blank(), 
          strip.text = element_blank())
})

# combine to one figure
gg_null <- cowplot::plot_grid(plotlist = gg_null_list, ncol = 2) |> 
  cowplot::ggdraw() +
  cowplot::draw_label("Stability value", x = 0.015, y = 0.5, angle = 90, size = base_size)

#### Create pop_n results ####

results_pop_df <- dplyr::left_join(x = results_null, y = mortality_df, 
                                     by = c("row_id", "pop_n", "nutrient_input")) |> 
  dplyr::mutate(include = dplyr::case_when(died_total > threshold_mort ~ "no", TRUE ~ "yes")) |> 
  dplyr::filter(include == "yes") |> 
  dplyr::mutate(part = factor(part, levels = c("ag_production", "bg_production", "ttl_production"),
                              labels = c("ag_production" = "Aboveground", "bg_production" = "Belowground", 
                                         "ttl_production" = "Total")), 
                measure = factor(measure, levels = c("alpha", "gamma", "beta", "synchrony"), 
                                 labels = c("Local", "Meta-Ecosys.", "Portfolio effect", "Synchrony"))) |> 
  dplyr::group_by(part, measure, nutrient_input, pop_n) |> #pop_n
  dplyr::summarise(mn = mean(value.cv), sd = sd(value.cv)) |> 
  dplyr::mutate(facet = dplyr::case_when(measure %in% c("Local", "Meta-Ecosys.") ~ "left", 
                                         TRUE ~ "right"))

gg_pop <- dplyr::filter(results_pop_df, part %in% c("Aboveground", "Total"), measure == "Portfolio effect") |> 
  ggplot() + 
  
  # adding geoms
  geom_hline(yintercept = 1, linetype = 2, color = 'grey') +
  geom_point(aes(x = pop_n, y = mn, color = nutrient_input)) + 
  geom_line(aes(x = pop_n, y = mn, color = nutrient_input, group = nutrient_input)) +
  geom_errorbar(aes(x = pop_n, ymin = mn - sd, ymax = mn + sd, color = nutrient_input), width = 0) +
  
  
  # adding facet labels
  geom_text(data = data.frame(x = c(0.75, 0.75), y = c(3.25, 2.75), part = c("Aboveground", "Total"),
                              label = c("i)", "ii)")), 
            aes(x = x , y = y, label = label), size = 5.5) +
  
  # facet wrap
  facet_wrap(. ~ part) + 
  
  # scales and labs
  scale_color_manual(name = "", values = color_enrich) +
  labs(x = "Population size", y = "Portfolio effect (PE)") +
  
  # adding box
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
  
  # facet_wrap by part
  facet_wrap(. ~ part, scales = "free") +
  
  # themeing
  theme_classic(base_size = base_size) +
  theme(legend.position = c(0.4, 0.925), legend.background = element_blank(),
        legend.text = element_text(size = base_size * 0.85),
        axis.line = element_blank(), strip.background = element_blank(), 
        strip.text = element_blank())
  
#### Create final plot

gg_combined <- cowplot::plot_grid(gg_null, gg_pop, nrow = 2,
                                  labels = c("A)", "B)"), label_fontface = "plain")

#### Save results #### 

suppoRt::save_ggplot(plot = gg_combined, filename = "Figure-S12.png",
                     path = "04_Figures/Supplemental/", width = width, height = height * 3 / 4,
                     units = units, dpi = dpi, overwrite = FALSE)

suppoRt::save_ggplot(plot = gg_combined, filename = "Figure-S12.pdf",
                     path = "04_Figures/Supplemental/", width = width, height = height * 3 / 4,
                     units = units, dpi = dpi, overwrite = FALSE)

