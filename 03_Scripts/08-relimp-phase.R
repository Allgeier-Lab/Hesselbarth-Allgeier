##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Variation partitioning phase variability

#### Load setup ####

source("05_Various/setup.R")

extension <- ".pdf"

#### Load/wrangle simulated data ####

amplitude <- "095"

file_path <- paste0("02_Data/05-variability-phase-", amplitude, ".rds")

df_results <- readr::read_rds(file_path) %>% 
  purrr::map_dfr(function(j) {
    dplyr::left_join(x = j$cv, y = j$production, by = c("part", "pop_n", "move_meta_sd", "phase_sd"), 
                     suffix = c(".cv", ".prod")) %>% 
      dplyr::mutate(value.move = mean(j$moved$moved, na.rm = TRUE))}) %>% 
  dplyr::mutate(part = factor(part, levels = c("ag_biomass", "bg_biomass", "ttl_biomass",
                                               "ag_production", "bg_production", "ttl_production")), 
                measure = factor(measure, levels = c("alpha", "gamma", "beta", "synchrony")),
                pop_n = factor(as.numeric(pop_n), ordered = TRUE)) %>% 
  dplyr::filter(part %in% c("ag_production", "bg_production", "ttl_production"), 
                measure %in% c("alpha", "gamma")) %>% 
  tibble::tibble()

#### Partition the variation ####

df_importance <- dplyr::group_by(df_results, part, measure, pop_n) %>% 
  dplyr::group_split() %>% 
  purrr::map_dfr(function(df_temp) {
    
    df_temp$value.cv <- log(df_temp$value.cv) - mean(log(df_temp$value.cv))
    
    lm_temp <- lm(value.cv ~ log(move_meta_sd) + log(phase_sd), data = df_temp)
    
    rel_r2 <- relaimpo::boot.relimp(lm_temp, type = "lmg", b = 1000, level = 0.95, fixed = FALSE) %>% 
      relaimpo::booteval.relimp(bty = "bca")
    
    tibble::tibble(
      part = unique(df_temp$part), measure = unique(df_temp$measure), pop_n = unique(df_temp$pop_n),
      beta = factor(x = c("connectivity", "phase", "residual"), levels  = c("residual", "phase", "connectivity")),
      mean = c(rel_r2@lmg, 1 - sum(rel_r2@lmg)), lower = c(rel_r2@lmg.lower, NA), higher = c(rel_r2@lmg.upper, NA)
    )
  })

#### setup ggplot globals ###

size_base <- 8.0

size_text <- 2.75

color_beta <- c(connectivity = "#f6c2a9", phase = "#b99dd9", residual = "grey")

color_part <- c(ag_production = "#df4e25", bg_production = "#007aa1", ttl_production = "#41b282")

#### Create ggplot ####

gg_indiv <- purrr::map(c(8, 16, 32, 64), function(pop_i) {
  
  purrr::map(c("alpha", "gamma"), function(measure_i) {
    
    df_importance_temp <- dplyr::filter(df_importance, pop_n == pop_i, measure == measure_i)
    
    gg_relimp_temp <- ggplot(data = df_importance_temp) + 
      geom_col(aes(x = part, y = mean, fill = beta)) + 
      scale_fill_manual(name = "", values = color_beta) +
      scale_color_manual(name = "", values = color_beta) +
      scale_x_discrete(labels = c(ag_production = "AG", bg_production = "BG", ttl_production = "Total"), 
                       limits = c("ttl_production", "bg_production", "ag_production")) +
      scale_y_continuous(limits = c(0, 1.0), n.breaks = 5) +
      labs(x = "", y = "Relative importance") +
      coord_flip() +
      theme_classic(base_size = size_base) + 
      theme(legend.position = "none", plot.margin = unit(c(t = 0.25, r = 0.25, b = 0.0, l = 0.25), "cm"))
    
    df_results_temp <- dplyr::filter(df_results, pop_n == pop_i, measure == measure_i) %>% 
      dplyr::select(part, move_meta_sd, phase_sd, value.cv) %>% 
      tidyr::pivot_longer(-c(part, value.cv))
    
    gg_lm <- ggplot(data = df_results_temp, aes(x = log(value), y = log(value.cv), color = part)) + 
      geom_point(shape = 1, alpha = 1) + 
      geom_smooth(size = 0.5, formula = y ~ x, se = FALSE, method = "lm") +
      scale_color_manual(name = "", values = color_part) +
      scale_y_continuous(breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 4), labels = function(x) round(x, 2)) +
      facet_wrap(. ~ name, labeller = labeller(name = c(move_meta_sd = "Connectivity", phase_sd = "Phase"))) + 
      labs(x = "log(Variability)", y = "log(CV)") +
      theme_classic(base_size = size_base) + 
      theme(legend.position = "none", strip.text = element_text(hjust = 0), strip.background = element_blank(),
            plot.margin = unit(c(t = 0.0, r = 0.25, b = ifelse(test = pop_i == 64, yes = 0.25, no = 0.75), l = 0.25), "cm"))
    
    cowplot::plot_grid(gg_relimp_temp, gg_lm, nrow = 2, rel_heights = c(0.3, 0.7))
    
  })}) %>% 
  purrr::flatten()

gg_dummy <- ggplot() + 
  geom_col(data = df_importance, aes(x = part, y = mean, fill = beta)) + 
  geom_line(data = df_results, aes(x = move_meta_sd, y = log(value.cv), color = part), size = 1.5) + 
  scale_fill_manual(name = "", values = color_beta, 
                    labels = c(connectivity = "Connectivity", phase = "Phase", residual = "Residuals")) +
  scale_color_manual(name = "", values = color_part, 
                     labels = c(ag_production = "Aboveground (AG)", bg_production = "Belowground (BG)", ttl_production = "Total")) +
  guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2)) +
  theme_classic(base_size = size_base) + 
  theme(legend.position = "bottom")

gg_varpart <- cowplot::plot_grid(cowplot::plot_grid(plotlist = gg_indiv, ncol = 2, nrow = 4, 
                                                    labels = "auto", label_fontface = "italic"), 
                                 cowplot::get_legend(gg_dummy), nrow = 2, rel_heights = c(0.975, 0.025))

suppoRt::save_ggplot(plot = gg_varpart, filename = paste0("08-varpart-phase-", amplitude, extension),
                     path = "04_Figures/", width = width, height = height,
                     units = units, dpi = dpi, overwrite = overwrite)

# ##--------------------------------------------##
# ##    Author: Maximilian H.K. Hesselbarth     ##
# ##    Coastal Ecology and Conservation Lab    ##
# ##    University of Michigan                  ##
# ##    mhessel@umich.edu                       ##
# ##    www.github.com/mhesselbarth             ##
# ##--------------------------------------------##
# 
# # Purpose: Variation partitioning phase variability
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
# df_results <- readr::read_rds(file_path) %>% 
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
# varpart_df <- dplyr::group_by(df_results, part, measure, pop_n) %>% 
#   dplyr::group_split() %>% 
#   purrr::map_dfr(function(df_temp) {
#     
#     part_temp <- vegan::varpart(log(df_temp$value.cv), ~ move_meta_sd, ~ phase_sd, data = df_temp)
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
