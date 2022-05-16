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

# nothing to change here

#### Stable values #### 

stable_values <- arrR::get_req_nutrients(bg_biomass = starting_list$bg_biomass,
                                         ag_biomass = starting_list$ag_biomass,
                                         parameters = parameters_list)

starting_list$nutrients_pool <- stable_values$nutrients_pool

starting_list$detritus_pool <- stable_values$detritus_pool

#### Setup experiment ####

experiment_df <- expand.grid(move_meta_mean = seq(from = 0, to = 1, by = 0.1),
                             move_meta_sd = seq(from = 0, to = 1, by = 0.1), 
                             amplitude_mn = amplitude_levels) %>% 
  dplyr::slice(rep(x = 1:dplyr::n(), each = 10)) %>% 
  tibble::tibble()

# table(experiment_df$move_meta_mean, experiment_df$move_meta_sd)

# setup metaecosystems
metasyst_temp <- meta.arrR::setup_meta(n = n, max_i = max_i, reef = reef_matrix,
                                       starting_values = starting_list, parameters = parameters_list,
                                       dimensions = dimensions, grain = grain, use_log = use_log, 
                                       verbose = FALSE)

foo_hpc <- function(move_meta_mean, move_meta_sd, amplitude_mn) {
  
  library(dplyr)
  
  # update move meta sd
  parameters_list$move_meta_mean <- move_meta_mean
  
  parameters_list$move_meta_sd <- move_meta_sd

  # create new attributed matrix
  attr_replace <- meta.arrR:::setup_attributes(fishpop = metasyst_temp$fishpop, parameters = parameters_list, 
                                               max_i = max_i)
  
  # replace matrix
  metasyst_temp$fishpop_attr[, 3] <- attr_replace[, 3]
  
  # simulate nutrient input
  input_temp <-  meta.arrR::simulate_nutr_input(n = n, max_i = max_i, 
                                                frequency = years, input_mn = nutrient_input, 
                                                amplitude_mn = amplitude_mn, amplitude_sd = 0.0,
                                                phase_mn = 0.0, phase_sd = 0.0,
                                                verbose = FALSE)
  
  # run model
  result_temp <- meta.arrR::run_simulation_meta(metasyst = metasyst_temp, parameters = parameters_list,
                                                nutrients_input = input_temp, movement = "behav",
                                                max_i = max_i, min_per_i = min_per_i,
                                                seagrass_each = seagrass_each,
                                                save_each = save_each, verbose = T)
  
  # filter model run
  result_temp <- meta.arrR::filter_meta(x = result_temp, filter = c((max_i / years) * years_filter, max_i), 
                                        reset = TRUE, verbose = FALSE)
  
  # get moved counts
  moved <- dplyr::bind_rows(result_temp$fishpop) %>% 
    dplyr::filter(timestep == max_i) %>% 
    dplyr::left_join(y = as.data.frame(metasyst_temp$fishpop_attr), by = "id") %>% 
    dplyr::select(id, moved, move_prob) %>% 
    dplyr::arrange(id)
  
  # calc cv
  cv <- meta.arrR::calc_variability(x = result_temp, lag = c(FALSE, TRUE))
  
  # calculate biomass/production
  production <- meta.arrR::summarize_meta(result = result_temp, biomass = TRUE, production = TRUE,
                                          lag = c(FALSE, FALSE)) %>% 
    purrr::map(tidyr::pivot_longer, cols = -c(meta, timestep), names_to = "part") %>% 
    dplyr::bind_rows() %>% 
    dplyr::filter(timestep == max_i) %>%
    dplyr::group_by(part) %>% 
    dplyr::summarise(alpha = mean(value), gamma = sum(value)) %>% 
    tidyr::pivot_longer(-part, names_to = "measure")
  
  # combine to result data.frame and list
  list(cv = dplyr::mutate(dplyr::bind_rows(cv), move_meta_mean = move_meta_mean, 
                          move_meta_sd = move_meta_sd, amplitude_mn = amplitude_mn), 
       production = dplyr::mutate(production, move_meta_mean = move_meta_mean, 
                                  move_meta_sd = move_meta_sd, amplitude_mn = amplitude_mn), 
       moved = dplyr::mutate(moved, move_meta_mean = move_meta_mean, 
                             move_meta_sd = move_meta_sd, amplitude_mn = amplitude_mn))
  
}

#### Submit HPC

globals <- c("parameters_list", "metasyst_temp", "max_i", "n", "nutrient_input",
             "min_per_i", "seagrass_each", "save_each", "years", "years_filter") 

sbatch_cv <- rslurm::slurm_apply(f = foo_hpc, params = experiment_df, 
                                 global_objects = globals, jobname = "cv_move_var",
                                 nodes = nrow(experiment_df), cpus_per_node = 1, 
                                 slurm_options = list("account" = account, 
                                                      "partition" = "standard",
                                                      "time" = "01:00:00", ## hh:mm::ss
                                                      "mem-per-cpu" = "7G", 
                                                      "exclude" = exclude_nodes),
                                 pkgs = c("dplyr", "meta.arrR"),
                                 rscript_path = rscript_path, sh_template = sh_template, 
                                 submit = FALSE)

#### Collect results #### 

suppoRt::rslurm_missing(x = sbatch_cv)

cv_result <- rslurm::get_slurm_out(sbatch_cv, outtype = "raw")

suppoRt::save_rds(object = cv_result, filename = "01-move-variability.rds",
                  path = "02_Data/", overwrite = overwrite)

rslurm::cleanup_files(sbatch_cv)

#### Load simulated data ####

move_var_list <- readr::read_rds("02_Data/01-move-variability.rds")

#### Setup ggplot ####

parts <- c("ag_production", "bg_production", "ttl_production")

measures <- c("alpha", "gamma")

titles_measure <- c(expression(paste("Local scale (", alpha, ")")),
                    expression(paste("Metaecosystem scale (", gamma, ")")))

titles_part <- c("Ag production", "Bg production", "Total production")

#### Create results figure CV ####

cv_result <- purrr::map_dfr(move_var_list, function(i) i$cv) %>% 
  dplyr::mutate(part = factor(part, levels = c("ag_biomass", "bg_biomass", "ttl_biomass", 
                                               "ag_production", "bg_production", "ttl_production")), 
                measure = factor(measure, levels = c("alpha", "gamma", "beta", "synchrony")), 
                amplitude_mn = factor(amplitude_mn, ordered = TRUE)) %>% 
  tibble::tibble()

for (i in amplitude_levels) {
  
  gg_cv <- purrr::map(seq_along(parts), function(j) {
    
    gg_measure <- purrr::map(seq_along(measures), function(k) {
      
      xlab_temp <- ifelse(test = j == 3, yes = expression(paste("Parameter ", italic(move_meta_mean))), 
                          no = "")
      
      ylab_temp <- ifelse(test = j == 2,
                          yes = expression(paste("Parameter ", italic(move_meta_sd))), 
                          no = "")
      
      xlab_temp_scat <- ifelse(test = j == 3, yes = "Parameter value", no = "")
      
      ylab_temp_scat <- ifelse(test = j == 2, yes = "Coeffiecent of variation", no = "")
      
      digits_temp <- ifelse(test = j == 1, yes = 2, no = ifelse(test = j == 2, yes = 3, no = 2))
      
      df_temp <- dplyr::filter(cv_result, part == parts[j], measure == measures[k],
                               amplitude_mn == i) %>% 
        dplyr::group_by(move_meta_mean, move_meta_sd) %>% 
        dplyr::summarise(value_mn = mean(value), value_sd = sd(value), .groups = "drop") %>% 
        dplyr::mutate(
          value_mn = dplyr::case_when(move_meta_mean == 0 & move_meta_sd == 0 ~ as.numeric(NA),
                                      TRUE ~ value_mn),
          value_sd = dplyr::case_when(move_meta_mean == 0 & move_meta_sd == 0 ~ as.numeric(NA),
                                      TRUE ~ value_sd),
          value_label = dplyr::case_when(value_mn == min(value_mn, na.rm = TRUE) ~ value_mn, 
                                         value_mn == max(value_mn, na.rm = TRUE) ~ value_mn, 
                                         TRUE ~ as.numeric(NA)), 
          label_color = dplyr::case_when(value_mn == min(value_mn, na.rm = TRUE) ~ "min", 
                                         value_mn == max(value_mn, na.rm = TRUE) ~ "max", 
                                         TRUE ~ as.character(NA))
          )
      
      row_median <- which(df_temp$value_mn == sort(df_temp$value_mn)[60])
      
      df_temp[row_median, "value_label"] <- df_temp[row_median, "value_mn"]
      
      df_temp[row_median, "label_color"] <-  "median"
      
      label_a = max(df_temp$value_mn, na.rm = TRUE) -
        (max(df_temp$value_mn, na.rm = TRUE) - min(df_temp$value_mn, na.rm = TRUE)) * 0.05

      label_b = max(df_temp$value_mn, na.rm = TRUE) -
        (max(df_temp$value_mn, na.rm = TRUE) - min(df_temp$value_mn, na.rm = TRUE)) * 0.2
      
      reg <- dplyr::select(df_temp, move_meta_mean, move_meta_sd, value_mn) %>% 
        tidyr::pivot_longer(-value_mn, names_to = "param_name", values_to = "param_value") %>% 
        dplyr::filter(!is.na(value_mn)) %>% 
        dplyr::group_by(param_name) %>% 
        dplyr::group_split() %>% 
        purrr::map(function(l) {
          
          lm_sum <- summary(lm(value_mn ~ param_value, data = l))
          
          c(s = signif(x = lm_sum$coefficients[[2]], digits = 1), r = signif(x = lm_sum$r.squared, digits = 2))
          
        })
      
      gg_scatter <- dplyr::select(df_temp, move_meta_mean, move_meta_sd, value_mn) %>% 
        tidyr::pivot_longer(-value_mn, names_to = "param_name", values_to = "param_value") %>% 
        dplyr::filter(!is.na(value_mn)) %>% 
        ggplot(aes(x = param_value, y = value_mn, color = param_name)) +
        geom_point(shape = 19, alpha = 0.25) +
        geom_smooth(formula = "y ~ x", method = "lm", se = FALSE) +
        annotate(geom = "label", x = 0.65, y = label_a, color = "#96c683", size = 2.5, parse = TRUE,
                 label = paste0("italic(mean):", "~slope==", reg[[1]][[1]], 
                                "~~R^{2}==", reg[[1]][[2]])) +
        annotate(geom = "label", x = 0.65, y = label_b, color = "#7f90e3", size = 2.5, parse = TRUE,
                 label = paste0("italic(sd):", "~slope==", reg[[2]][[1]], 
                                "~~R^{2}==", reg[[2]][[2]])) +
        scale_color_manual(name = "", values = c("#96c683", "#7f90e3")) + 
        labs(x = xlab_temp_scat, y = ylab_temp_scat) +
        theme_classic(base_size = 10) +
        theme(legend.position = "none", 
              plot.margin = unit(c(0.25, 0.5, 0.25, 0.25), "cm")) # t, r, b, l
    
      gg_tile <- ggplot(data = df_temp, aes(x = move_meta_mean, y = move_meta_sd, fill = value_mn)) + 
        geom_tile() + 
        geom_label(aes(label = round(value_label, digits_temp), color = label_color), size = 3) +
        labs(x = xlab_temp, y = ylab_temp) +
        scale_fill_viridis_c(option = "A", na.value = "grey") + 
        scale_color_manual(values = c("black", "grey", "white")) + 
        scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
        scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
        coord_equal() +
        theme_classic(base_size = 10) + 
        theme(legend.position = "none",
              axis.line = element_blank(), panel.border = element_rect(size = 1, fill = NA), 
              aspect.ratio = 1, plot.margin = unit(c(0.25, 0.25, 0.5, 0.25), "cm")) # t, r, l, b
      
      cowplot::plot_grid(gg_tile, gg_scatter, ncol = 2, rel_widths = c(0.5, 0.5))
      
    })
    
    cowplot::plot_grid(plotlist = gg_measure, nrow = 1, ncol = 2)
    
  })
  
  gg_cv <- cowplot::plot_grid(plotlist = gg_cv, nrow = 3, ncol = 1, labels = c("a)", "b)", "c)"))
  
  ii <- ifelse(i == 0.05, yes = "low", no = ifelse(i == 0.5, yes = "med", no = "hi"))

  suppoRt::save_ggplot(plot = gg_cv, filename = paste0("01-move-variability-cv-", ii, ".pdf"), 
                       path = "04_Figures/", width = height, height = width, units = units,  dpi = dpi, 
                       overwrite = overwrite)
}

#### Create results figure production ####

prod_result <- purrr::map_dfr(move_var_list, function(i) {
  
  prod_temp <- dplyr::filter(i$production, part %in% c("ag_production", "bg_production", "ttl_production"))
  
  cv_temp <- dplyr::filter(i$cv, part %in% c("ag_production", "bg_production", "ttl_production"), 
                measure %in% c("alpha", "gamma"))
  
  dplyr::left_join(x = prod_temp, y = cv_temp, by = c("part", "measure", "move_meta_mean", 
                                                      "move_meta_sd", "amplitude_mn"), 
                   suffix = c(".prd", ".cv")) %>% 
    dplyr::mutate(value.mvd = mean(i$moved$moved))}) %>%
  dplyr::mutate(part = factor(part, levels = c("ag_biomass", "bg_biomass", "ttl_biomass", 
                                               "ag_production", "bg_production", "ttl_production")), 
                measure = factor(measure, levels = c("alpha", "gamma", "beta", "synchrony")), 
                amplitude_mn = factor(amplitude_mn, ordered = TRUE)) %>% 
  dplyr::relocate(value.prd, .before = value.cv) %>% 
  tibble::tibble()

for (i in amplitude_levels) {
  
  gg_prod <- purrr::map(seq_along(parts), function(j) {
    
    gg_measure <- purrr::map(seq_along(measures), function(k) {
      
      xlab_temp_move <- ifelse(test = j == 3, yes = "Mean cross-system movement", no = "")
      
      xlab_temp_cv <- ifelse(test = j == 3, yes = "Coeffiecent of variation", no = "")
      
      ylab_temp <- ifelse(test = j == 2, yes = "Total primary production", no = "")
      
      df_temp <- dplyr::filter(prod_result, part == parts[j], measure == measures[k], 
                               amplitude_mn == i, value.mvd != 0.0)
      
      ylabel_annotate <- max(df_temp$value.prd, na.rm = TRUE) -
        (max(df_temp$value.prd, na.rm = TRUE) - min(df_temp$value.prd, na.rm = TRUE)) * 0.05
      
      xlabel_annotate <- max(df_temp$value.cv, na.rm = TRUE) -
        (max(df_temp$value.cv, na.rm = TRUE) - min(df_temp$value.cv, na.rm = TRUE)) * 0.25
      
      reg_mod <- summary(lm(value.prd ~ value.cv, data = df_temp))
      
      ggplot(data = df_temp, aes(x = value.cv, y = value.prd)) + 
        geom_point(shape = 1, alpha = 0.25) + 
        geom_smooth(se = FALSE, linetype = 1, color = "#007aa1", 
                    formula = y ~ x, method = "lm") + 
        annotate(geom = "label", x = xlabel_annotate, y = ylabel_annotate,
                 color = "#007aa1", size = 2.5, parse = TRUE,
                 label = paste0("~slope==", signif(reg_mod$coefficients[[2]], digits = 1),
                                "~~R^{2}==", signif(reg_mod$r.squared, digits = 2))) +
        labs(x = xlab_temp_cv, y = "") +
        theme_classic(base_size = 10) + 
        theme(plot.margin = unit(c(0.25, 0.25, 0.5, 0.25), "cm"))
      
    })
    
    cowplot::plot_grid(plotlist = gg_measure, nrow = 1, ncol = 2)
    
  })
  
  gg_prod <- cowplot::plot_grid(plotlist = gg_prod, nrow = 3, ncol = 1, 
                                labels = c("a)", "b)", "c)"))
  
  ii <- ifelse(i == 0.05, yes = "low", no = ifelse(i == 0.5, yes = "med", no = "hi"))
  
  suppoRt::save_ggplot(plot = gg_prod, filename = paste0("01-move-variability-prod-", ii, ".pdf"), 
                       path = "04_Figures/", width = width, height = height, units = units,  dpi = dpi, 
                       overwrite = overwrite)
}

#### Create results figure comparison production ####

prod_result_comp <- purrr::map_dfr(move_var_list, function(i) i$production) %>% 
  dplyr::mutate(part = factor(part, levels = c("ag_biomass", "bg_biomass", "ttl_biomass", 
                                               "ag_production", "bg_production", "ttl_production")), 
                measure = factor(measure, levels = c("alpha", "gamma", "beta", "synchrony")), 
                amplitude_mn = factor(amplitude_mn, ordered = TRUE, labels = c("Low", "Medium", "High"))) %>% 
  tibble::tibble()

gg_prod_comp <- purrr::map(seq_along(parts), function(j) {
    
  gg_measure <- purrr::map(seq_along(measures), function(k) {
    
    legend_temp <- ifelse(test = c(j == 1, j == 1), yes = c(0.85, 0.2), no = c(0.15, 0.2))
    
    xlab_temp <- ifelse(test = j == 3, yes = "Amplitude treatment", no = "")
    
    ylab_temp <- ifelse(test = k == 1, yes = "Total primary production", no = "")
    
    guides_temp <- ifelse(test = j == 3, yes = "legend", no = "none")
    
    df_temp <- dplyr::filter(prod_result_comp, part == parts[j], measure == measures[k], 
                             move_meta_mean != 0 & move_meta_sd != 0)
    
    ggplot(data = df_temp, aes(x = amplitude_mn, y = value, 
                               color = amplitude_mn, fill = amplitude_mn)) + 
      geom_violin(alpha = 0.5) +
      geom_boxplot(width = 0.25, outlier.shape = NA, fill = NA) +
      scale_color_manual(name = "Amplitude", values = c("#b63116", "#fad248", "#339cca")) + 
      scale_fill_manual(name = "Amplitude", values = c("#b63116", "#fad248", "#339cca")) +
      guides(color = "none", fill = guides_temp) +
      labs(x = xlab_temp, y = ylab_temp) + 
      theme_classic(base_size = 10) + 
      theme(legend.position = legend_temp, legend.text = element_text(size = 7.5),
            legend.title = element_text(size = 7.5), 
            legend.box.background = element_rect(colour = "black"))
    
  })
  
  cowplot::plot_grid(plotlist = gg_measure, nrow = 1, ncol = 2)
    
})
  
gg_prod_comp <- cowplot::plot_grid(plotlist = gg_prod_comp, nrow = 3, ncol = 1, 
                                   labels = c("a)", "b)", "c)"))
  
suppoRt::save_ggplot(plot = gg_prod_comp, filename = paste0("01-move-variability-prod-comp.pdf"), 
                     path = "04_Figures/", width = width, height = height, units = units,  dpi = dpi, 
                     overwrite = overwrite)
