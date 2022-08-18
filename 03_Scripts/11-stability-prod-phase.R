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
source("05_Various/import_data.R")

#### Load/wrangle simulated data ####

amplitude <- "095"

df_results <- import_data(path = paste0("02_Data/05-variability-phase-", amplitude, ".rds"))

#### Fit regression model ####

df_regression <- dplyr::filter(df_results, part %in% c("ag_production", "bg_production", "ttl_production"), 
                               measure %in% c("alpha", "gamma")) %>%
  dplyr::group_by(part, measure, pop_n) %>%
  dplyr::group_split() %>%
  purrr::map_dfr(function(df_temp) {
    
    df_temp$value.prod <- log(df_temp$value.prod) - mean(log(df_temp$value.prod))
    
    lm_temp <- lm(value.prod ~ log(value.cv), data = df_temp)
    
    broom::tidy(lm_temp) %>%
      dplyr::mutate(part = unique(df_temp$part), measure = unique(df_temp$measure),
                    pop_n = unique(df_temp$pop_n), .before = term) %>% 
      dplyr::mutate(term = dplyr::case_when(term == "(Intercept)" ~ "Intercept", TRUE ~ term), 
                    p.value = dplyr::case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**",
                                               p.value < 0.05 ~  "*", p.value >= 0.05 ~ "n.s."),
                    direction = dplyr::case_when(term != "Intercept" & estimate < 0.0 ~ "increase", 
                                                 term != "Intercept" & estimate > 0.0 ~ "decrease", 
                                                 term == "Intercept" ~ as.character(NA)))
    
  })

#### Setup ggplot ####

colors_pop <- c("8" = "#447861", "16" = "#13315f", "32" = "#59386c", "64" = "#b24422")

size_point <- 1.0
size_line <- 0.75
size_base <- 10.0

alpha <- 0.25

#### Create ggplot ####

gg_scatter_list <- dplyr::filter(df_results, measure %in% c("alpha", "gamma")) %>%
  dplyr::group_by(part, measure) %>%
  dplyr::group_split() %>%
  purrr::map(function(df_temp) {
    
    # init ggplot
    ggplot(data = df_temp, aes(x = log(value.cv), y = log(value.prod), color = pop_n)) + 
      
      # adding geoms
      geom_point(shape = 1, alpha = alpha, size = size_point) + 
      geom_smooth(formula = y ~ x, method = "lm", se = FALSE, size = size_line) +
      
      # facet wrap 
      facet_wrap(. ~ pop_n, scales = "free", nrow = 2, ncol = 2) + 
      
      # set scales
      scale_color_manual(name = "Population size", values = colors_pop) +
      scale_x_continuous(limits = function(x) range(x), labels = function(x) round(exp(x), 2),
                         breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 4)) +
      scale_y_continuous(limits = function(x) range(x), labels = function(x) round(exp(x), 1),
                         breaks = function(x) seq(quantile(x, 0.1), quantile(x, 0.9), length.out = 4)) +
      
      # labels and themes
      labs(x = "", y = "") +
      theme_classic(base_size = 10.0) + 
      theme(legend.position = "none", plot.title = element_text(size = 8.0), 
            strip.background = element_blank(), strip.text = element_blank(),
            plot.margin = margin(t = 5.5, r = 5.5, b = 0.0, l = 5.5, unit = "pt"))
    
  })

gg_dummy <- ggplot(data = dplyr::filter(df_results, measure != "beta"), 
                   aes(x = value.cv, y = value.prod, color = pop_n)) + 
  geom_point(shape = 1, alpha = alpha) + 
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE, size = size_line) +
  scale_color_manual(name = "Population size", values = colors_pop) +
  theme_classic(base_size = size_base) + 
  theme(legend.position = "bottom")

gg_combined <- cowplot::plot_grid(plotlist = gg_scatter_list, nrow = 3, ncol = 2, byrow = TRUE, 
                                  labels = c("a)", "b)", "c)", "d)", "e)", "f)"),
                                  label_fontface = "plain", label_size = 12.0)

gg_combined <- cowplot::ggdraw(gg_combined, xlim = c(-0.015, 1.0)) +
  cowplot::draw_label("Coeffiecent of variation", x = 0.5, y = 0, vjust = -0.5, angle = 0, size = size_base) + 
  cowplot::draw_label("Primary production", x = 0.0, y = 0.5, vjust = 0.0, angle = 90, size = size_base)

gg_combined <- cowplot::plot_grid(gg_combined, cowplot::get_legend(gg_dummy),
                                  nrow = 2, ncol = 1, rel_heights = c(0.95, 0.05))

#### Save plot ####

suppoRt::save_ggplot(plot = gg_combined, filename = paste0("11-phase-", amplitude, extension),
                     path = "04_Figures/", width = height, height = width * 0.75,
                     units = units, dpi = dpi, overwrite = FALSE)
