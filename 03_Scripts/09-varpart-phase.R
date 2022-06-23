##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Variation partitioning

#### Load setup ####

source("05_Various/setup.R")

extension <- ".pdf"

#### Load/wrangle simulated data ####

amplitude <- "095"

file_path <- paste0("02_Data/05-variability-phase-", amplitude, ".rds")

df_cv_prod <- readr::read_rds(file_path) %>% 
  purrr::map_dfr(function(j) {
    dplyr::left_join(x = j$cv, y = j$production, by = c("part", "pop_n", "move_meta_sd", "phase_sd"), 
                     suffix = c(".cv", ".prod")) %>% 
      dplyr::mutate(value.move = mean(j$moved$moved, na.rm = TRUE))
  }) %>% 
  dplyr::mutate(part = factor(part, levels = c("ag_biomass", "bg_biomass", "ttl_biomass",
                                               "ag_production", "bg_production", "ttl_production")), 
                measure = factor(measure, levels = c("alpha", "gamma", "beta", "synchrony")),
                pop_n = factor(as.numeric(pop_n), ordered = TRUE)) %>% 
  dplyr::filter(part %in% c("ag_production", "bg_production", "ttl_production"), 
                measure %in% c("alpha", "gamma", "beta")) %>% 
  tibble::tibble()

#### Partition the variation ####

varpart_df <- dplyr::group_by(df_cv_prod, part, measure, pop_n) %>% 
  dplyr::group_split() %>% 
  purrr::map_dfr(function(df_temp) {
    
    part_temp <- vegan::varpart(log(df_temp$value.cv), ~ move_meta_sd, ~ phase_sd, data = df_temp)
    
    tibble::tibble(part = unique(df_temp$part), measure = unique(df_temp$measure),
                   pop_n = unique(df_temp$pop_n),
                   explanatory = c("connectivity", "shared", "phase", "residuals"),
                   x = c(0.5, 2.0, 3.5, 5.0), y = c(1, 1, 1, 1),
                   value = part_temp$part$indfract$Adj.R.square) %>% 
      dplyr::mutate(value = dplyr::case_when(value < 0 ~ 0, TRUE ~ value))}) %>% 
  dplyr::mutate(y = dplyr::case_when(part == "ag_production" ~ y + 0.45,
                                     part == "bg_production" ~ y, 
                                     part == "ttl_production" ~ y - 0.45))

df_circle <- tibble::tibble(x = c(1, 3), y = c(1, 1), 
                            explanatory = c("connectivity", "phase"))

#### setup ggplot globals ###

base_size <- 8.0

text_size <- 2.75

color_var <- c(connectivity = "#f6c2a9", phase = "#b99dd9")

color_part <- c(ag_production = "#df4e25", bg_production = "#007aa1", ttl_production = "#41b282")

#### Create ggplot ####

gg_indiv <- purrr::map(c(8, 16, 32, 64), function(pop_i) {
  
  gg_measure <- purrr::map(c("alpha", "gamma"), function(measure_i) {
    
    varpart_df_temp <- dplyr::filter(varpart_df, pop_n == pop_i, measure == measure_i)
    
    df_cv_prod_temp <- dplyr::filter(df_cv_prod, pop_n == pop_i, measure == measure_i) %>% 
      dplyr::select(part, move_meta_sd, phase_sd, value.cv) %>% 
      tidyr::pivot_longer(-c(part, value.cv))
    
    # strip_text <- element_blank()
    
    # if (pop_i == 8) strip_text <- element_text(hjust = 0)
    
    gg_varpar_temp <- ggplot(data = varpart_df_temp) + 
      ggforce::geom_circle(data = df_circle, aes(x0 = x, y0 = y, r = 1.5, fill = explanatory),
                           alpha = 0.75, color = "black", size = 0.25) +
      annotate(geom = "text", x = 0.5, y = 2, label = "Connectivity", color = "black", size = text_size) +
      annotate(geom = "text", x = 3.5, y = 2, label = "Phase", color = "black", size = text_size) +
      annotate(geom = "text", x = 5.0, y = 2, label = "Resid.", size = text_size) +
      geom_text(aes(x = x, y = y, color = part, label = round(value, 2)), size = text_size) +
      scale_color_manual(name = "", values = color_part) +
      scale_fill_manual(values = color_var) +
      guides(fill = "none") +
      coord_fixed(ratio = 0.5) +
      theme_void(base_size = base_size) +
      theme(legend.position = "none")
    
    gg_lm <- ggplot(data = df_cv_prod_temp, aes(x = value, y = log(value.cv), color = part)) + 
      geom_point(shape = 1, alpha = 0.15) + 
      geom_smooth(size = 0.5, formula = y ~ x, se = FALSE, method = "lm") +
      scale_color_manual(name = "", values = color_part) +
      scale_y_continuous(breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 4), 
                         labels = function(x) round(x, 2)) +
      facet_wrap(. ~ name, labeller = labeller(name = c(move_meta_sd = "Connectivity",
                                                        phase_sd = "Phase"))) + 
      labs(x = ifelse(test = pop_i == 64, yes = "Variability", no = ""), 
           y = ifelse(test = measure_i == "alpha", yes = "log(CV)", no = "")) +
      theme_classic(base_size = base_size) + 
      theme(legend.position = "none", 
            strip.text = element_blank(), strip.background = element_blank(), 
            axis.title = element_text(size = 6.5))
    
    cowplot::plot_grid(gg_varpar_temp, gg_lm, nrow = 2, rel_heights = c(0.35, 0.65))
    
  })
  
  cowplot::plot_grid(plotlist = gg_measure, ncol = 2)
  
})

gg_dummy <- ggplot(data = df_cv_prod, aes(x = move_meta_sd, y = value.cv, color = part)) +
  geom_point() + 
  geom_smooth(formula = y ~ x, se = FALSE, method = "lm") + 
  scale_color_manual(name = "", values = color_part, labels = c("Aboveground", "Belowground", "Total")) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom", 
        strip.text = element_text(hjust = 0), strip.background = element_blank())

gg_varpart <- cowplot::plot_grid(cowplot::plot_grid(plotlist = gg_indiv, nrow = 4), 
                                 cowplot::get_legend(gg_dummy), nrow = 2, rel_heights = c(0.975, 0.025))

suppoRt::save_ggplot(plot = gg_varpart, filename = paste0("09-varpart-phase-", amplitude, extension),
                     path = "04_Figures/", width = width, height = height,
                     units = units, dpi = dpi, overwrite = T)

# ##--------------------------------------------##
# ##    Author: Maximilian H.K. Hesselbarth     ##
# ##    Coastal Ecology and Conservation Lab    ##
# ##    University of Michigan                  ##
# ##    mhessel@umich.edu                       ##
# ##    www.github.com/mhesselbarth             ##
# ##--------------------------------------------##
# 
# # Purpose: Variation partitioning
# 
# #### Load setup ####
# 
# source("05_Various/setup.R")
# 
# extension <- ".pdf"
# 
# #### Load/wrangle simulated data ####
# 
# amplitude <- "095"
# 
# file_path <- paste0("02_Data/05-variability-phase-", amplitude, ".rds")
# 
# df_cv_prod <- readr::read_rds(file_path) %>% 
#   purrr::map_dfr(function(j) {
#     dplyr::left_join(x = j$cv, y = j$production, by = c("part", "pop_n", "move_meta_sd", "phase_sd"), 
#                      suffix = c(".cv", ".prod")) %>% 
#       dplyr::mutate(value.move = mean(j$moved$moved, na.rm = TRUE))
#   }) %>% 
#   dplyr::mutate(part = factor(part, levels = c("ag_biomass", "bg_biomass", "ttl_biomass",
#                                                "ag_production", "bg_production", "ttl_production")), 
#                 measure = factor(measure, levels = c("alpha", "gamma", "beta", "synchrony")),
#                 pop_n = factor(as.numeric(pop_n), ordered = TRUE)) %>% 
#   dplyr::filter(part %in% c("ag_production", "bg_production", "ttl_production"), 
#                 measure %in% c("alpha", "gamma", "beta")) %>% 
#   tibble::tibble()
# 
# #### Partition the variation ####
# 
# varpart_df <- dplyr::group_by(df_cv_prod, part, measure, pop_n) %>% 
#   dplyr::group_split() %>% 
#   purrr::map_dfr(function(df_temp) {
#     
#     part_temp <- vegan::varpart(df_temp$value.cv, ~ move_meta_sd, ~ phase_sd, data = df_temp)
#     
#     tibble::tibble(part = unique(df_temp$part), measure = unique(df_temp$measure),
#                    pop_n = unique(df_temp$pop_n),
#                    explanatory = c("connectivity", "shared", "phase", "residuals"),
#                    x = c(0.5, 2.0, 3.5, 4.5), y = c(1, 1, 1, -0.5),
#                    value = part_temp$part$indfract$Adj.R.square) %>% 
#       dplyr::mutate(value = dplyr::case_when(value < 0 ~ 0, TRUE ~ value))}) %>% 
#   dplyr::mutate(y = dplyr::case_when(part == "ag_production" ~ y + 0.25,
#                                      part == "bg_production" ~ y, 
#                                      part == "ttl_production" ~ y - 0.25))
# 
# df_circle <- tibble::tibble(x = c(1, 3), y = c(1, 1), 
#                             explanatory = c("connectivity", "phase"))
# 
# #### setup ggplot globals ###
# 
# text_size <- 3.5
# 
# color_var <- c(connectivity = "#f6c2a9", phase = "#b99dd9")
# 
# color_part <- c(ag_production = "#df4e25", bg_production = "#007aa1", ttl_production = "#41b282")
# 
# #### Create ggplot ####
# 
# gg_varpart <- dplyr::filter(varpart_df, measure != "beta") %>% 
#   ggplot() + 
#   ggforce::geom_circle(data = df_circle, aes(x0 = x, y0 = y, r = 1.5, fill = explanatory),
#                        alpha = 0.75, color = "black", size = 0.25) +
#   annotate(geom = "text", x = 0.5, y = 1.75, label = "Connectivity", color = "black", size = text_size) +
#   annotate(geom = "text", x = 3.5, y = 1.75, label = "Phase", color = "black", size = text_size) +
#   annotate(geom = "text", x = 4.0, y = -0.5, label = "Resd.", size = text_size) +
#   geom_text(aes(x = x, y = y, color = part, label = round(value, 2)), size = text_size) +
#   scale_color_manual(name = "", values = color_part) +
#   scale_fill_manual(values = color_var) +
#   facet_grid(rows = dplyr::vars(pop_n), cols = dplyr::vars(measure), 
#              labeller = labeller(pop_n = c("8" = "8 Indiv.", "16" = "16 Indiv.", 
#                                            "32" = "32 Indiv.", "64" = "64 Indiv."), 
#                                  measure = c(alpha = "Local scale", gamma = "Meta-ecosystem scale"))) + 
#   guides(fill = "none") +
#   coord_equal() +
#   theme_void() + theme(legend.position = "bottom")
# 
# suppoRt::save_ggplot(plot = gg_varpart, filename = paste0("09-varpart-phase-", amplitude, extension),
#                      path = "04_Figures/", width = width, height = height,
#                      units = units, dpi = dpi, overwrite = overwrite)
