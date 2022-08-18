##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Create total result figure including regression parameters and relative importance

#### Load setup ####

source("05_Various/setup.R")
source("05_Various/import_data.R")

#### Load/wrangle simulated data ####

amplitude <- "095"

df_results <- import_data(path = paste0("02_Data/05-variability-phase-", amplitude, ".rds"))

#### Split data #### 

df_list <- dplyr::filter(df_results, part %in% c("ag_production", "bg_production", "ttl_production"), 
              measure %in% c("alpha", "gamma")) %>%
  dplyr::group_by(part, measure, pop_n) %>%
  dplyr::group_split()

#### Fit regression model ####

df_regression <- purrr::map_dfr(df_list, function(df_temp) {

    df_temp$value.cv <- log(df_temp$value.cv) - mean(log(df_temp$value.cv))

    lm_temp <- lm(value.cv ~ log(biotic) * log(abiotic), data = df_temp)

    broom::tidy(lm_temp) %>%
      dplyr::mutate(part = unique(df_temp$part), measure = unique(df_temp$measure),
                    pop_n = unique(df_temp$pop_n), .before = term) %>% 
      dplyr::mutate(term = c("Intercept", "Biotic", "Abiotic", "Interaction"))
  })

#### Relative importance R2 ####

df_importance <- purrr::map_dfr(df_list, function(df_temp) {
  
  df_temp$value.cv <- log(df_temp$value.cv) - mean(log(df_temp$value.cv))
  
  lm_temp <- lm(value.cv ~ log(biotic) * log(abiotic), data = df_temp)
  
  rel_r2 <- relaimpo::boot.relimp(lm_temp, type = "lmg", b = 1000, level = 0.95, fixed = FALSE) %>%
    relaimpo::booteval.relimp(bty = "bca")
  
  tibble::tibble(
    part = unique(df_temp$part), measure = unique(df_temp$measure), pop_n = unique(df_temp$pop_n),
    beta = factor(x = c("Biotic", "Abiotic", "Interaction", "Residual"),
                  levels  = c("Residual", "Biotic", "Abiotic", "Interaction")),
    mean = c(rel_r2@lmg, 1 - sum(rel_r2@lmg)), lower = c(rel_r2@lmg.lower, NA), higher = c(rel_r2@lmg.upper, NA),
    p_value = ifelse(c(summary(lm_temp)$coefficients[c(2:4), 4], NA) < 0.05, yes = "signif.", no = "n.s.")
  )
})

#### Setup ggplot ####

size_base <- 10.0
size_text <- 2.0
size_point <- 3.5

color_parameter <- c("Intercept" = "#df4e25", "Biotic" = "#41b282", "Abiotic" = "#007aa1", "Interaction" = "#fcb252")
color_relimp <- c("Biotic" = "#41b282", "Abiotic" = "#007aa1", "Interaction" = "#fcb252", "Residual" = "grey")

#### Create ggplot ####

gg_list <- purrr::map(c("alpha", "gamma"), function(measure_i) {
  
  purrr::map(c("ag_production", "bg_production", "ttl_production"), function(part_i) {
    
    df_regression_temp <- dplyr::filter(df_regression, measure == measure_i, part == part_i) %>% 
      dplyr::select(-c(std.error, statistic)) %>% 
      dplyr::mutate(p.value = dplyr::case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**",
                                               p.value < 0.05 ~  "*", p.value >= 0.05 ~ ""))
    
    df_importance_temp <- dplyr::filter(df_importance, measure == measure_i, part == part_i)
    
    w <- 0.5
    
    gg_regression <- ggplot(data = df_regression_temp) + 
      
      # zero line
      geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
      
      # Lines
      geom_line(aes(x = pop_n, y = estimate, group = term, color = term),
                alpha = 0.5, position = position_dodge(width = w)) +
      
      # Points
      geom_point(aes(x = pop_n, y = estimate, color = term),
                 size = size_point, shape = 19, position = position_dodge(width = w)) +
      
      # Text
      geom_text(aes(x = pop_n, y = estimate, label = p.value, group = term), color = "white",
                size = size_text, position = position_dodge(width = w), vjust = 0.75) +
      
      # set scales
      scale_color_manual(name = "Scale", values = color_parameter) +
      scale_y_continuous(limits = c(min(df_regression$estimate), max(df_regression$estimate)),
                         breaks = seq(min(df_regression$estimate), max(df_regression$estimate), length.out = 4), 
                         labels = function(x) round(x, digits = 2)) + 
      coord_cartesian(clip = "off") +
      
      # labels and themes
      labs(x = "", y = "") +
      theme_classic(base_size = 10.0) + 
      theme(legend.position = "none", plot.title = element_text(size = 8.0), 
            strip.background = element_blank(), strip.text = element_blank(), 
            plot.margin = margin(t = 5.5, r = -5.5, b = 5.5, l = 5.5, unit = "pt"))
    
    gg_relimp <- ggplot(data = df_importance_temp) + 
      
      # relative importance bars
      geom_col(aes(x = pop_n, y = mean * 100, fill = beta)) + 
      
      # set scales
      scale_fill_manual(name = "", values = color_relimp) +
      scale_y_continuous(labels = function(x) paste0(x, "%")) + 
      
      # labels and themes
      labs(x = "", y = "") +
      theme_classic(base_size = size_base) + 
      theme(legend.position = "none", plot.title = element_text(size = 8.0), 
            strip.background = element_blank(), strip.text = element_blank(), 
            plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = -5.5, unit = "pt"))
    
    cowplot::plot_grid(gg_regression, gg_relimp, ncol = 2, rel_widths = c(0.5, 0.5))
    
  })
}) %>% purrr::flatten()

#### Combine plots ####

gg_dummy <- data.frame(beta = factor(c("Biotic", "Abiotic", "Interaction", "Intercept", "Residual"), 
                                     levels = c("Biotic", "Abiotic", "Interaction", "Intercept", "Residual")),
                       mean = c(1, 1, 1, 1, 1)) %>% 
  ggplot() + 
  geom_col(aes(x = beta, y = mean, fill = beta)) + 
  scale_fill_manual(name = "", values = c(color_relimp["Biotic"], color_relimp["Abiotic"], color_relimp["Interaction"], 
                                          color_parameter["Intercept"], color_relimp["Residual"])) +
  guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2)) +
  theme_classic(base_size = size_base) + 
  theme(legend.position = "bottom")

gg_combined <- cowplot::plot_grid(plotlist = gg_list, nrow = 3, ncol = 2, byrow = FALSE, 
                                  labels = c("a)", "b)", "c)", "d)", "e)", "f)"), 
                                  label_fontface = "plain", label_size = 12.0)

gg_combined <- cowplot::ggdraw(gg_combined, xlim = c(-0.015, 1.0)) + 
  cowplot::draw_label("Population size", x = 0.5, y = 0, vjust = -0.5, angle = 0, size = size_base) + 
  cowplot::draw_label("Parameter estimate / Relative importance [%]", x = 0.0, y = 0.5, vjust = 0.0, 
                      angle = 90, size = size_base)

gg_combined <- cowplot::plot_grid(gg_combined, cowplot::get_legend(gg_dummy),
                                  nrow = 2, ncol = 1, rel_heights = c(0.95, 0.05))

#### Save plot ####

suppoRt::save_ggplot(plot = gg_combined, filename = paste0("08-phase-", amplitude, extension),
                     path = "04_Figures/", width = height, height = width * 0.7,
                     units = units, dpi = dpi, overwrite = FALSE)
