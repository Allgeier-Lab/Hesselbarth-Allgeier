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

#### Relative importance R2 ####

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
      mean = c(rel_r2@lmg, 1 - sum(rel_r2@lmg)), lower = c(rel_r2@lmg.lower, NA), higher = c(rel_r2@lmg.upper, NA), 
      p_value = ifelse(c(summary(lm_temp)$coefficients[c(2:3), 4], NA) < 0.05, yes = "signif.", no = "n.s.")
    )
  })

#### setup ggplot globals ###

size_base <- 10.0

size_text <- 2.75

color_beta <- c(connectivity = "#88a0dc", phase = "#f9d14a", residual = "grey")

color_part <- c(ag_production = "#df4e25", bg_production = "#007aa1", ttl_production = "#41b282")

cv_minmax <- df_results$value.cv %>% 
  log() %>% 
  range()

#### Create relative importance ggplot ####

gg_relimp_list <- purrr::map(c("alpha", "gamma"), function(measure_i) {
  
  df_importance_temp <- dplyr::filter(df_importance, measure == measure_i)
  
  ggplot(data = df_importance_temp) + 
    geom_col(aes(x = pop_n, y = mean, fill = beta)) + 
    facet_wrap(. ~ part, nrow = 3,scales = "free", 
               labeller = labeller(part = c("ag_production" = ifelse(measure_i == "alpha", yes = "Aboveground", no = ""), 
                                            "bg_production" = ifelse(measure_i == "alpha", yes = "Belowground", no = ""), 
                                            "ttl_production" = ifelse(measure_i == "alpha", yes = "Total", no = "")))) +
    scale_fill_manual(name = "", values = color_beta) +
    scale_y_continuous(limits = c(0, 1.0), n.breaks = 5) +
    labs(title = ifelse(measure_i == "alpha", yes = "Local scale", no = "Meta-ecosystem scale"), 
         x = "Population size", y = ifelse(measure_i == "alpha", yes = "Relative importance", no = "")) +
    theme_classic(base_size = size_base) + 
    theme(legend.position = "none", plot.title = element_text(size = 8.0), 
          strip.background = element_blank(), strip.text = element_text(hjust = 0))
})

gg_relimp_dummy <- ggplot() + 
  geom_col(data = df_importance, aes(x = part, y = mean, fill = beta)) + 
  scale_fill_manual(name = "", values = color_beta, 
                    labels = c(connectivity = "Biotic variability", phase = "Abiotic variability", residual = "Residuals")) +
  guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2)) +
  theme_classic(base_size = size_base) + 
  theme(legend.position = "bottom")

gg_relimp_ttl <- cowplot::plot_grid(cowplot::plot_grid(plotlist = gg_relimp_list, ncol = 2), 
                                    cowplot::get_legend(gg_relimp_dummy), nrow = 2, rel_heights = c(0.975, 0.025))

suppoRt::save_ggplot(plot = gg_relimp_ttl, filename = paste0("08-relimp-phase-", amplitude, extension),
                     path = "04_Figures/", width = width, height = height,
                     units = units, dpi = dpi, overwrite = overwrite)
