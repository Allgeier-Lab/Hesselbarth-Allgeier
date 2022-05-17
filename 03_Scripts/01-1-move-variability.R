##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

#### Load setup ####

source("05_Various/setup.R")
source("01_Functions/convert_notification.R")

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
                                                save_each = save_each, verbose = FALSE)
  
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

parts <- c("ag_production", "bg_production", "ttl_production")

measures <- c("alpha", "gamma")

#### Create results figure CV ####

cv_result <- purrr::map_dfr(move_var_list, function(i) {
    dplyr::mutate(i$cv, moved_mn = mean(i$moved$moved, na.rm = TRUE), .before = "amplitude_mn")}) %>% 
  dplyr::mutate(part = factor(part, levels = c("ag_biomass", "bg_biomass", "ttl_biomass", 
                                               "ag_production", "bg_production", "ttl_production")), 
                measure = factor(measure, levels = c("alpha", "gamma", "beta", "synchrony")), 
                amplitude_mn = factor(amplitude_mn, ordered = TRUE)) %>% 
  tibble::tibble()

for (i in amplitude_levels) {
  
  gg_cv <- purrr::map(seq_along(parts), function(j) {
    
    gg_measure <- purrr::map(seq_along(measures), function(k) {
      
      # setup labels #
      
      xlab_tile <- ifelse(test = j == 3, yes = expression(paste("Parameter ", italic(move_meta_mean))), no = "")
      
      ylab_tile <- ifelse(test = j == 2, yes = expression(paste("Parameter ", italic(move_meta_sd))), no = "")
      
      xlab_scatter <- ifelse(test = j == 3, yes = "Mean cross-ecosystem movement", no = "")
      
      ylab_scatter <- ifelse(test = j == 2, yes = "Coeffiecent of variation", no = "")
      
      digits_tile <- ifelse(test = j == 1, yes = 2, no = ifelse(test = j == 2, yes = 3, no = 2))
      
      # data tile plot #
      
      df_tile <- dplyr::filter(cv_result, part == parts[[j]], measure == measures[[k]],
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
                                         TRUE ~ as.character(NA)))
      
      row_median <- which(df_tile$value_mn == sort(df_tile$value_mn)[60])
      
      df_tile[row_median, "value_label"] <- df_tile[row_median, "value_mn"]
      
      df_tile[row_median, "label_color"] <-  "median"
      
      # ggplot tile #
      
      gg_tile <- ggplot(data = df_tile, aes(x = move_meta_mean, y = move_meta_sd, fill = value_mn)) + 
        geom_tile() + 
        geom_label(aes(label = round(value_label, digits_tile), color = label_color), size = 2.0) +
        labs(x = xlab_tile, y = ylab_tile) +
        scale_fill_viridis_c(option = "A", na.value = "grey") + 
        scale_color_manual(values = c("black", "grey", "white")) + 
        scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
        scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
        coord_equal() +
        theme_classic(base_size = 8) + 
        theme(legend.position = "none",
              axis.line = element_blank(), panel.border = element_rect(size = 1, fill = NA), 
              aspect.ratio = 1, plot.margin = unit(c(0.25, 0.25, 0.5, 0.25), "cm")) # t, r, l, b

      # data scatter plot #    
    
      df_scatter <- dplyr::filter(cv_result, part == parts[[j]], measure == measures[[k]],
                                  amplitude_mn == i, moved_mn != 0)
      
      x_text <- min(df_scatter$moved_mn) + (max(df_scatter$moved_mn) - min(df_scatter$moved_mn)) * 0.5
      
      y_text <- min(df_scatter$value) + (max(df_scatter$value) - min(df_scatter$value)) * 0.95
      
      # regression model #
      
      lm_sum <- summary(lm(value ~ moved_mn, data = df_scatter))
      
      lm_intercept <- ifelse(lm_sum$coefficients[[1]] < 0.001 || lm_sum$coefficients[[1]] > 1000, 
                             yes = convert_notification(lm_sum$coefficients[[1]]), 
                             no = round(lm_sum$coefficients[[1]], digits = 3))
      
      lm_coef <- ifelse(lm_sum$coefficients[[2]] < 0.001 || lm_sum$coefficients[[2]] > 1000, 
                        yes = convert_notification(abs(lm_sum$coefficients[[2]])), 
                             no = round(abs(lm_sum$coefficients[[2]]), digits = 3))
      
      lm_rr <- ifelse(lm_sum$r.squared < 0.001 || lm_sum$r.squared > 1000, 
                      yes = convert_notification(lm_sum$r.squared), 
                      no = round(lm_sum$r.squared, digits = 3))
      
      lm_direction <- ifelse(lm_sum$coefficients[[2]] > 0, yes = "+", no = "-")
    
      ## ggplot scatter plot #
      
      gg_scatter <- ggplot(data = df_scatter, aes(x = moved_mn, y = value)) +
        geom_point(shape = 19, alpha = 0.15) +
        geom_smooth(formula = "y ~ x", method = "lm", se = FALSE, color = "#007aa1") + 
        annotate(geom = "label", x = x_text, y = y_text, size = 2.0, color = "#007aa1", parse = TRUE,
                 label = paste0("italic(y)==", lm_intercept, lm_direction, lm_coef, "*italic(x)", 
                                "*';'~~italic(R)^{2}==", lm_rr)) +
        scale_x_continuous(breaks = seq(from = min(df_scatter$moved_mn), to = max(df_scatter$moved_mn), 
                                        length.out = 5)) +
        labs(x = xlab_scatter, y = ylab_scatter) +
        theme_classic(base_size = 8) +
        theme(legend.position = "none", 
              plot.margin = unit(c(0.25, 0.5, 0.25, 0.25), "cm")) # t, r, b, l
      
      # combine ggplots #
      
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
    dplyr::mutate(value.mvd = mean(i$moved$moved, na.rm = TRUE), .before = "amplitude_mn")}) %>%
  dplyr::mutate(part = factor(part, levels = c("ag_biomass", "bg_biomass", "ttl_biomass", 
                                               "ag_production", "bg_production", "ttl_production")), 
                measure = factor(measure, levels = c("alpha", "gamma", "beta", "synchrony")), 
                amplitude_mn = factor(amplitude_mn, ordered = TRUE)) %>% 
  dplyr::relocate(value.prd, .before = value.cv) %>% 
  tibble::tibble()

for (i in amplitude_levels) {
  
  gg_prod <- purrr::map(seq_along(parts), function(j) {
    
    gg_measure <- purrr::map(seq_along(measures), function(k) {
      
      # setup labels #
      
      xlab_cv <- ifelse(test = j == 3, yes = "Coeffiecent of variation", no = "")
      
      ylab_cv <- ifelse(test = j == 2, yes = "Total primary production", no = "")

      # data scatter plot #
      
      df_scatter <- dplyr::filter(prod_result, part == parts[j], measure == measures[k], 
                                  amplitude_mn == i, value.mvd != 0.0)
      
      x_text <- min(df_scatter$value.cv) + (max(df_scatter$value.cv) - min(df_scatter$value.cv)) * 0.5
      
      y_text <- min(df_scatter$value.prd) + (max(df_scatter$value.prd) - min(df_scatter$value.prd)) * 0.95
      
      # regression model #
      
      lm_sum <- summary(lm(value.prd ~ value.cv, data = df_scatter))
      
      lm_intercept <- ifelse(lm_sum$coefficients[[1]] < 0.001 || lm_sum$coefficients[[1]] > 1000, 
                             yes = convert_notification(lm_sum$coefficients[[1]]), 
                             no = round(lm_sum$coefficients[[1]], digits = 3))
      
      lm_coef <- ifelse(lm_sum$coefficients[[2]] < 0.001 || lm_sum$coefficients[[2]] > 1000, 
                        yes = convert_notification(abs(lm_sum$coefficients[[2]])), 
                        no = round(abs(lm_sum$coefficients[[2]]), digits = 3))
      
      lm_rr <- ifelse(lm_sum$r.squared < 0.001 || lm_sum$r.squared > 1000, 
                      yes = convert_notification(lm_sum$r.squared), 
                      no = round(lm_sum$r.squared, digits = 3))
      
      lm_direction <- ifelse(lm_sum$coefficients[[2]] > 0, yes = "+", no = "-")
      
      ## ggplot scatter plot #
      
      ggplot(data = df_scatter, aes(x = value.cv, y = value.prd)) + 
        geom_point(shape = 1, alpha = 0.25) + 
        geom_smooth(se = FALSE, linetype = 1, color = "#007aa1", 
                    formula = y ~ x, method = "lm") + 
        annotate(geom = "label", x = x_text, y = y_text, size = 2.0, color = "#007aa1", parse = TRUE,
                 label = paste0("italic(y)==", lm_intercept, lm_direction, lm_coef, "*italic(x)",
                                "*';'~~italic(R)^{2}==", lm_rr)) +
        labs(x = xlab_cv, y = ylab_cv) +
        theme_classic(base_size = 8) + 
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
    
    # setup labels #
    
    legend_comp <- ifelse(test = c(j == 1, j == 1), yes = c(0.85, 0.2), no = c(0.15, 0.2))
    
    xlab_comp <- ifelse(test = j == 3, yes = "Amplitude treatment", no = "")
    
    ylab_comp <- ifelse(test = k == 1, yes = "Total primary production", no = "")
    
    guides_comp <- ifelse(test = j == 3, yes = "legend", no = "none")
    
    # get scatterplot #
    
    df_comp <- dplyr::filter(prod_result_comp, part == parts[j], measure == measures[k], 
                             move_meta_mean != 0 & move_meta_sd != 0)
    
    # ggplot scatterplot #
    
    ggplot(data = df_comp, aes(x = amplitude_mn, y = value, 
                               color = amplitude_mn, fill = amplitude_mn)) + 
      geom_violin(alpha = 0.5) +
      geom_boxplot(width = 0.25, outlier.shape = NA, fill = NA) +
      scale_color_manual(name = "Amplitude", values = c("#b63116", "#fad248", "#339cca")) + 
      scale_fill_manual(name = "Amplitude", values = c("#b63116", "#fad248", "#339cca")) +
      guides(color = "none", fill = guides_comp) +
      labs(x = xlab_comp, y = ylab_comp) + 
      theme_classic(base_size = 10) + 
      theme(legend.position = legend_comp, legend.text = element_text(size = 7.5),
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

#### Create PE  ####

pe_result <- purrr::map_dfr(move_var_list, function(i) {
    dplyr::mutate(i$cv, moved = mean(i$moved$moved, na.rm = TRUE), .before = "amplitude_mn")}) %>%
  dplyr::mutate(part = factor(part, levels = c("ag_biomass", "bg_biomass", "ttl_biomass",
                                               "ag_production", "bg_production", "ttl_production")),
                measure = factor(measure, levels = c("alpha", "gamma", "beta", "synchrony")),
                amplitude_mn = factor(amplitude_mn, ordered = TRUE, labels = c("Low", "Medium", "High"))) %>%
  tibble::tibble()

gg_pe <- purrr::map(seq_along(parts), function(j) {
  
  # setup labels #
  
  xlab_pe <- ifelse(test = j == 3, yes = "Mean cross-ecosystem movement", no = "")
  
  ylab_pe <- ifelse(test = j == 2, yes = "Portfolio effect", no = "")

  # get data scatter plot #
  
  df_pe <- dplyr::filter(pe_result, part == parts[j], measure == "beta",
                           move_meta_mean != 0 & move_meta_sd != 0)
  
  range_pe <- dplyr::filter(pe_result, measure == "beta") %>%
    dplyr::pull(value) %>%
    range()
  
  label_x <- min(df_pe$moved, na.rm = TRUE) + (max(df_pe$moved, na.rm = TRUE) - min(df_pe$moved, na.rm = TRUE)) * 0.85
  
  label_y_lo <- range_pe[[1]] + (range_pe[[2]] - range_pe[[1]]) * 0.95
  
  label_y_med <- range_pe[[1]] + (range_pe[[2]] - range_pe[[1]]) * 0.85
  
  label_y_hi <- range_pe[[1]] + (range_pe[[2]] - range_pe[[1]]) * 0.75
  
  # regression model #
  
  lm_sum <- dplyr::group_by(df_pe, amplitude_mn) %>% 
    dplyr::group_split() %>% 
    purrr::map(function(l) {
      lm_sum <- summary(lm(value ~ moved, data = l))
      
      lm_intercept <- ifelse(lm_sum$coefficients[[1]] < 0.001 || lm_sum$coefficients[[1]] > 1000, 
                             yes = convert_notification(lm_sum$coefficients[[1]]), 
                             no = round(lm_sum$coefficients[[1]], digits = 3))
      
      lm_coef <- ifelse(lm_sum$coefficients[[2]] < 0.001 || lm_sum$coefficients[[2]] > 1000, 
                        yes = convert_notification(abs(lm_sum$coefficients[[2]])), 
                        no = round(abs(lm_sum$coefficients[[2]]), digits = 3))
      
      lm_rr <- ifelse(lm_sum$r.squared < 0.001 || lm_sum$r.squared > 1000, 
                      yes = convert_notification(lm_sum$r.squared), 
                      no = round(lm_sum$r.squared, digits = 3))
      
      lm_direction <- ifelse(lm_sum$coefficients[[2]] > 0, yes = "+", no = "-")
      
      c(int = lm_intercept, coef = lm_coef, rr = lm_rr, direction = lm_direction)})
  
  # gg scatterplot #
  
  ggplot(data = df_pe, aes(x = moved, y = value, color = amplitude_mn)) + 
    geom_point(alpha = 0.25) + 
    geom_smooth(se = FALSE, linetype = 1, formula = y ~ x, method = "lm") + 
    geom_hline(yintercept = 1, linetype = 2, color =  "grey") + 
    annotate(geom = "label", x = label_x, y = label_y_lo, size = 2.0, color = "#b63116", parse = TRUE,
             label = paste0("italic(y)==", lm_sum[[1]]["int"], lm_sum[[1]]["direction"], 
                            lm_sum[[1]]["coef"], "*italic(x)", "*';'~~italic(R)^{2}==", lm_sum[[1]]["rr"])) +
    annotate(geom = "label", x = label_x, y = label_y_med, size = 2.0, color = "#fad248", parse = TRUE,
             label = paste0("italic(y)==", lm_sum[[2]]["int"], lm_sum[[2]]["direction"], 
                            lm_sum[[2]]["coef"], "*italic(x)", "*';'~~italic(R)^{2}==", lm_sum[[2]]["rr"])) +
    annotate(geom = "label", x = label_x, y = label_y_hi, size = 2.0, color = "#339cca", parse = TRUE,
             label = paste0("italic(y)==", lm_sum[[3]]["int"], lm_sum[[3]]["direction"], 
                            lm_sum[[3]]["coef"], "*italic(x)", "*';'~~italic(R)^{2}==", lm_sum[[3]]["rr"])) +
    scale_y_continuous(breaks = seq(from = 0, to = ceiling(range_pe[[2]]), by = 1),
                       limits = c(0, ceiling(range_pe[[2]]))) +
    scale_x_continuous(breaks = seq(from = min(df_pe$moved), to = max(df_pe$moved), length.out = 5)) +
    scale_color_manual(name = "Amplitude", values = c("#b63116", "#fad248", "#339cca")) +
    guides(color = "none") +
    labs(x = xlab_pe, y = ylab_pe) +
    theme_classic(base_size = 8) + 
    theme(legend.position = c(0.15, 0.075), legend.background = element_blank())

})

gg_pe <- cowplot::plot_grid(plotlist = gg_pe, nrow = 3, ncol = 1, 
                            labels = c("a)", "b)", "c)"))

suppoRt::save_ggplot(plot = gg_pe, filename = paste0("01-move-variability-pe.pdf"), 
                     path = "04_Figures/", width = width, height = height, units = units,  dpi = dpi, 
                     overwrite = T)
