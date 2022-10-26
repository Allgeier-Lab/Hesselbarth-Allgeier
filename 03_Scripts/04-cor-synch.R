##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Calculate correlations between variability, cv, and synchrony

#### Load setup ####

source("05_Various/setup.R")
source("01_Functions/import_data.R")

#### Load/wrangle simulated data ####

n <- 5

results_phase_df <- import_data(path = paste0("02_Data/result-phase-", n, ".rds"))

results_noise_df <- import_data(path = paste0("02_Data/result-noise-", n, ".rds"))

results_combined_df <- dplyr::bind_rows(phase = results_phase_df, noise = results_noise_df, 
                                        .id = "scenario") |>
  dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise")))

#### ggplot globals ####

size_base <- 10

digits <- 2

color_scale_var <- c("biotic" = "#41b282", "abiotic" = "#007aa1")

# color_scale_scl <- c("alpha" = "#cbcbcb", "gamma" = "#32b2da")

gg_dummy_var <- data.frame(beta = c("biotic", "abiotic"),
                           mean = c(1, 1)) |> 
  dplyr::mutate(beta = factor(beta, levels = c("biotic", "abiotic"))) |> 
  ggplot() + 
  geom_point(aes(x = beta, y = mean, color = beta)) + 
  geom_line(aes(x = beta, y = mean, color = beta)) + 
  scale_color_manual(name = "", values = color_scale_var, 
                    labels = c("Consumer behavior", "Abiotic subsidies")) +
  guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2)) +
  theme_classic(base_size = size_base) + 
  theme(legend.position = "bottom")

# gg_dummy_scl <- data.frame(beta = c("alpha", "gamma"),
#                            mean = c(1, 1)) |> 
#   dplyr::mutate(beta = factor(beta, levels = c("alpha", "gamma"))) |> 
#   ggplot() + 
#   geom_point(aes(x = beta, y = mean, color = beta)) + 
#   geom_line(aes(x = beta, y = mean, color = beta)) + 
#   scale_color_manual(name = "", values = color_scale_scl, 
#                      labels = c("Local scale", "Meta-ecosystem scale")) +
#   guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2)) +
#   theme_classic(base_size = size_base) + 
#   theme(legend.position = "bottom")

#### Reshape data ####

factor_lvls <- paste(rep(x = c("ag", "bg", "ttl"), each = 5), rep(x = c(8, 16, 32, 64, 128), times = 3), sep = "_")

results_var_synch_df <- dplyr::select(results_combined_df, scenario, part, measure, pop_n, biotic, abiotic, value.cv) |> 
  dplyr::filter(part %in% c("ag_production", "bg_production", "ttl_production"), measure == "synchrony") |> 
  tidyr::pivot_longer(cols = c(biotic, abiotic), names_to = "variability", values_to = "value.var") |> 
  tidyr::unite(col = "treatment", part, pop_n, remove = FALSE) |> 
  dplyr::mutate(treatment = stringr::str_remove(treatment, pattern = "_production"), 
                treatment = factor(treatment, levels = factor_lvls))

# results_cv_synch_df <- dplyr::select(results_combined_df, scenario, part, measure, pop_n, value.cv) |> 
#   dplyr::filter(part %in% c("ag_production", "bg_production", "ttl_production"), measure != "beta") |> 
#   dplyr::mutate(row_id = rep(x = 1:(dplyr::n() / 3), each = 3)) |> 
#   tidyr::pivot_wider(names_from = measure, values_from = value.cv) |> 
#   tidyr::pivot_longer(cols = c(alpha, gamma), names_to = "measure") |> 
#   tidyr::unite(col = "treatment", part, pop_n, remove = FALSE) |> 
#   dplyr::mutate(treatment = stringr::str_remove(treatment, pattern = "_production"), 
#                 treatment = factor(treatment, levels = factor_lvls))

#### Create ggplot variability vs. synchrony ####

gg_var_synch <- purrr::map(c("phase", "noise"), function(i) {
  
  gg_part <- purrr::map(c("ag_production", "bg_production", "ttl_production"), function(j) {
    
    data_temp <- dplyr::filter(results_var_synch_df, scenario == i, part == j)
      
    if (j == "ag_production") label_part <-  "Aboveground"
    if (j == "bg_production") label_part <-  "Belowground"
    if (j == "ttl_production") label_part <- "Total      "
    
    p <- ggplot(data = data_temp, aes(x = value.var, y = value.cv, color = variability)) + 
      
      # adding geoms
      geom_point(shape = 1, alpha = 0.5) + 
      geom_smooth(method = "lm", formula = "y ~ x", se = FALSE, size = 0.5) +
      
      # adding labels
      ggpmisc::stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..adj.rr.label..,
                                                                    sep = "~~~")),
                            parse = TRUE, coef.digits = 2, rr.digits = 2, size = 2.0, 
                            geom = "label_npc", alpha = 1.0, vstep = 0.2, hstep = 0.0) +
      
      # facet wrap 
      facet_wrap(. ~ pop_n, ncol = 5, scales = "free_y", 
                 labeller = labeller(pop_n = function(x) paste("Pop. size", x))) +
      
      # scales
      scale_x_continuous(breaks = seq(from = 0.1, to = 1, by = 0.25)) +
      scale_y_continuous(breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 4),
                         labels = function(x) round(x, 2)) +
      scale_color_manual(name = "", values = color_scale_var) +
      
      # themes
      labs(x = "", y = "") +
      # labs(title = label_part, subtitle = label_pop_n) +
      theme_classic(base_size = size_base) + 
      theme(strip.background = element_blank(), strip.text = element_text(hjust = 0), 
            legend.position = "none",
            plot.margin = unit(c(t = 0.0, r = 0.0, b = 0.0, l = 0.0), "mm"))
    
    cowplot::ggdraw(p, xlim = c(0.0, 1.05)) + 
      cowplot::draw_label(label = label_part, x = 1.0, y = 0.65,
                          vjust = -0.5, angle = 270, size = size_base * 0.85)
      
  })
  
  gg_combined <- cowplot::plot_grid(plotlist = gg_part, nrow = 3)
  
  gg_combined <- cowplot::ggdraw(gg_combined, xlim = c(-0.05, 1.0), ylim = c(-0.05, 1.0)) + 
    cowplot::draw_label("Variability", x = 0.5, y = 0, vjust = 0.5, angle = 0, size = size_base) + 
    cowplot::draw_label("Synchrony", x = 0.0, y = 0.5, vjust = -0.5, angle = 90, size = size_base)
  
  cowplot::plot_grid(gg_combined, cowplot::get_legend(gg_dummy_var),
                     nrow = 2, ncol = 1, rel_heights = c(0.95, 0.05))
  
})

names(gg_var_synch) <- c("phase", "noise")

# ##### Create ggplot synchrony vs. cv ####
# 
# gg_synch_cv <- purrr::map(c("phase", "noise"), function(i) {
#   
#   gg_part <- purrr::map(c("ag_production", "bg_production", "ttl_production"), function(j) {
#     
#     data_temp <- dplyr::filter(results_cv_synch_df, scenario == i, part == j)
#       
#     if (j == "ag_production") label_part <-  "Aboveground"
#     if (j == "bg_production") label_part <-  "Belowground"
#     if (j == "ttl_production") label_part <- "Total      "
#     
#     p <- ggplot(data = data_temp, aes(x = synchrony, y = value, color = measure)) + 
#       
#       # adding geoms
#       geom_point(shape = 1, alpha = 0.15) + 
#       geom_smooth(method = "lm", formula = "y ~ x", se = FALSE, size = 0.5) +
#       
#       # adding labels
#       ggpmisc::stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..adj.rr.label..,
#                                                                sep = "~~~")),
#                             parse = TRUE, coef.digits = 2, rr.digits = 2, size = 2.25) +
#       
#       # facet wrap 
#       facet_wrap(. ~ pop_n, ncol = 5, scales = "free", 
#                  labeller = labeller(pop_n = function(x) paste("Pop. size", x))) +
#       
#       # scales
#       scale_x_continuous(breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 4),
#                          labels = function(x) round(x, 2)) +
#       scale_y_continuous(breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 4),
#                          labels = function(x) round(x, 2)) +
#       scale_color_manual(name = "", values = color_scale_scl, 
#                          labels = c("Local scale", "Meta-ecosystem scale")) +
#       
#       # themes
#       labs(x = "", y = "") +
#       theme_classic(base_size = 10) + 
#       theme(strip.background = element_blank(), strip.text = element_text(hjust = 0.0), 
#             legend.position = "none")
#     
#     cowplot::ggdraw(p, xlim = c(0.0, 1.05)) + 
#       cowplot::draw_label(label = label_part, x = 1.0, y = 0.65,
#                           vjust = -0.5, angle = 270, size = size_base * 0.85)
#     
#   })
#   
#   gg_combined <- cowplot::plot_grid(plotlist = gg_part, nrow = 3)
#   
#   gg_combined <- cowplot::ggdraw(gg_combined, xlim = c(-0.05, 1.0), ylim = c(-0.05, 1.0)) + 
#     cowplot::draw_label("Synchrony", x = 0.5, y = 0, vjust = 0.5, 
#                         angle = 0, size = size_base) + 
#     cowplot::draw_label("Coefficient of variation", x = 0.0, y = 0.5, vjust = -0.5,
#                         angle = 90, size = size_base)
#   
#   cowplot::plot_grid(gg_combined, cowplot::get_legend(gg_dummy_scl),
#                      nrow = 2, ncol = 1, rel_heights = c(0.95, 0.05))
#   
# })
# 
# names(gg_synch_cv) <- c("phase", "noise")

##### Save plots ####

suppoRt::save_ggplot(plot = gg_var_synch$noise, filename = paste0("Figure-2", extension),
                     path = "04_Figures/", width = height, height = width * 0.75,
                     units = units, dpi = dpi, overwrite = T)

# suppoRt::save_ggplot(plot = gg_var_synch$phase, filename = paste0("Figure-2", extension),
#                      path = "04_Figures/Appendix/", width = height, height = width * 0.75,
#                      units = units, dpi = dpi, overwrite = FALSE)
# 
# ##--------------------------------------------##
# ##    Author: Maximilian H.K. Hesselbarth     ##
# ##    Coastal Ecology and Conservation Lab    ##
# ##    University of Michigan                  ##
# ##    mhessel@umich.edu                       ##
# ##    www.github.com/mhesselbarth             ##
# ##--------------------------------------------##
# 
# # Purpose: Calculate correlations between variability, cv, and synchrony
# 
# #### Load setup ####
# 
# source("05_Various/setup.R")
# source("01_Functions/import_data.R")
# 
# #### Load/wrangle simulated data ####
# 
# n <- 5
# 
# results_phase_df <- import_data(path = paste0("02_Data/result-phase-", n, ".rds"))
# 
# results_noise_df <- import_data(path = paste0("02_Data/result-noise-", n, ".rds"))
# 
# results_combined_df <- dplyr::bind_rows(phase = results_phase_df, noise = results_noise_df, 
#                                         .id = "scenario") |>
#   dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise")))
# 
# #### ggplot globals ####
# 
# size_base <- 10
# 
# digits <- 2
# 
# color_scale_var <- c("biotic" = "#41b282", "abiotic" = "#007aa1")
# 
# color_scale_scl <- c("alpha" = "#cbcbcb", "gamma" = "#32b2da")
# 
# gg_dummy_var <- data.frame(biotic = seq(0, 1, 0.1), abiotic = seq(0, 1, 0.1), 
#                            value.cv = runif(n = 11)) |> 
#   ggplot() + 
#   geom_raster(aes(x = biotic, y = abiotic, fill = value.cv)) + 
#   scale_fill_gradientn(name = "", colors = rev(MetBrewer::met.brewer(name = "Hokusai1", n = 255)), 
#                        limits = c(0,1), breaks = c(0, 0.5, 1.0), labels = c("low", "med", "high")) +
#   theme_classic(base_size = size_base) + 
#   theme(legend.position = "bottom", 
#         legend.key.width = unit(15, "mm"), legend.key.height = unit(2.5, 'mm'))
# 
# #### Reshape data ####
# 
# factor_lvls <- paste(rep(x = c("ag", "bg", "ttl"), each = 5), rep(x = c(8, 16, 32, 64, 128), times = 3), 
#                      sep = "_")
# 
# results_var_synch_df <- dplyr::select(results_combined_df, scenario, part, measure, pop_n, biotic, abiotic, value.cv) |> 
#   dplyr::filter(part %in% c("ag_production", "bg_production", "ttl_production"), measure == "synchrony") |> 
#   tidyr::pivot_longer(cols = c(biotic, abiotic), names_to = "variability", values_to = "value.var") |> 
#   tidyr::unite(col = "treatment", part, pop_n, remove = FALSE) |> 
#   dplyr::mutate(treatment = stringr::str_remove(treatment, pattern = "_production"), 
#                 treatment = factor(treatment, levels = factor_lvls))
# 
# #### Create ggplot variability vs. synchrony ####
# 
# gg_var_synch <- purrr::map(c("phase", "noise"), function(i) {
#   
#   gg_part <- purrr::map(c("ag_production", "bg_production", "ttl_production"), function(j) {
#     
#     data_temp <- dplyr::filter(results_var_synch_df, scenario == i, part == j) |> 
#       tidyr::pivot_wider(names_from = variability, values_from = value.var) |> 
#       dplyr::mutate(biotic = cut(biotic, breaks = seq(0, 1, 0.2) , include.lowest = TRUE), 
#                     abiotic = cut(abiotic, breaks = seq(0, 1, 0.2), include.lowest = TRUE)) |> 
#       dplyr::group_by(biotic, abiotic, pop_n) |> 
#       dplyr::summarise(value.cv = mean(value.cv), .groups = "drop") |>
#       dplyr::group_by(pop_n) |> 
#       dplyr::mutate(value.cv = scales::rescale(value.cv), .groups = "drop")
#     
#     if (j == "ag_production") label_part <-  "Aboveground"
#     if (j == "bg_production") label_part <-  "Belowground"
#     if (j == "ttl_production") label_part <- "Total      "
#     
#     # if (j == "ag_production") label_pop_n <- paste0("Population size: ", k)
#     
#     p <- ggplot(data = data_temp, aes(x = biotic, y = abiotic, fill = value.cv)) + 
#       
#       # adding geoms
#       geom_raster() +
#       
#       # facet wrap
#       facet_wrap(. ~ pop_n, ncol = 5, labeller = labeller(pop_n = function(x) 
#         paste("Population size", x))) +
#       
#       # scales
#       labs(x = "", y = "") +
#       scale_x_discrete(labels = c("<0.2", "", "<0.6", "", "<1.0")) +
#       scale_y_discrete(labels = c("<0.2", "", "<0.6", "", "<1.0")) +
#       scale_fill_gradientn(colors = rev(MetBrewer::met.brewer(name = "Hokusai1", n = 255)), 
#                            limits = c(0, 1)) +
#       coord_fixed(ratio = 1) +
#       
#       # themes
#       # labs(title = label_part) +
#       theme_classic(base_size = size_base) +
#       theme(strip.background = element_blank(), strip.text = element_text(hjust = 0), 
#             legend.position = "none",
#             axis.line = element_blank(), panel.background = element_rect(fill = NA, color = "black"))
#     
#     cowplot::ggdraw(p, xlim = c(0.0, 1.05)) + 
#       cowplot::draw_label(label = label_part, x = 1.0, y = 0.65,
#                           vjust = -0.5, angle = 270, size = size_base * 0.85)
#     
#   })
#   
#   gg_combined <- cowplot::plot_grid(plotlist = gg_part, nrow = 3)
#   
#   gg_combined <- cowplot::ggdraw(gg_combined, xlim = c(-0.05, 1.0), ylim = c(-0.05, 1.0)) + 
#     cowplot::draw_label("Consumer behavior", x = 0.5, y = 0, vjust = 0.5, angle = 0, size = size_base) + 
#     cowplot::draw_label("Abiotic subsidies", x = 0.0, y = 0.5, vjust = -0.5, angle = 90, size = size_base)
#   
#   cowplot::plot_grid(gg_combined, cowplot::get_legend(gg_dummy_var),
#                      nrow = 2, ncol = 1, rel_heights = c(0.95, 0.05))
#   
# })
# 
# names(gg_var_synch) <- c("phase", "noise")
# 
# ##### Save plots ####
# 
# # extension <- ".png"
# 
# suppoRt::save_ggplot(plot = gg_var_synch$noise, filename = paste0("Figure-2-raster", extension),
#                      path = "04_Figures/", width = height, height = width * 0.75,
#                      units = units, dpi = dpi, overwrite = FALSE)
# 
# # suppoRt::save_ggplot(plot = gg_var_synch$phase, filename = paste0("Figure-2-alt", extension),
# #                      path = "04_Figures/Appendix/", width = height, height = width * 0.75,
# #                      units = units, dpi = dpi, overwrite = FALSE)

