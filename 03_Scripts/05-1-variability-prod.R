##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Create total result figure of variability and cumulative PP

#### Load setup ####

source("01_Functions/setup.R")
source("01_Functions/import-cv.R")
source("01_Functions/import-mortality.R")

#### Load/wrangle simulated data ####

near <- FALSE

results_phase_df <- import_cv(path = "02_Data/result-phase.rds", near = near)

results_noise_df <- import_cv(path = "02_Data/result-noise.rds", near = near)

results_combined_df <- dplyr::bind_rows(phase = results_phase_df, noise = results_noise_df, 
                                        .id = "scenario") |> 
  dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise")))

#### Load mortality data ####

mortality_phase_df <- import_mortality(path = "02_Data/result-phase.rds") |> 
  dplyr::filter(pop_n > 0)

mortality_noise_df <- import_mortality(path = "02_Data/result-noise.rds") |> 
  dplyr::filter(pop_n > 0)

mortality_combined_df <- dplyr::bind_rows(phase = mortality_phase_df, noise = mortality_noise_df,
                                          .id = "scenario") |>
  dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise"))) |>
  dplyr::group_by(scenario, row_id, pop_n, nutrient_input) |>
  dplyr::summarise(died_total = mean(died_total), .groups = "drop")

#### Filter abundances/mortality #### 

results_combined_df <- dplyr::left_join(x = results_combined_df, y = mortality_combined_df, 
                                        by = c("scenario", "row_id", "pop_n", "nutrient_input")) |> 
  dplyr::mutate(include = dplyr::case_when(pop_n == 128 & nutrient_input == "low" &
                                             died_total > threshold_mort ~ "no", 
                                           TRUE ~ "yes")) |>
  dplyr::filter(pop_n != 0, abiotic != 0.0, biotic != 0.0, include == "yes")

#### Split data #### 

results_combined_list <- dplyr::left_join(x = dplyr::filter(results_combined_df, measure == "beta") |> 
                                            dplyr::select(scenario, row_id, part, pop_n, nutrient_input, value.cv),
                                          y = dplyr::filter(results_combined_df, measure == "gamma") |> 
                                            dplyr::select(scenario, row_id, part, pop_n, nutrient_input, value.prod), 
                                          by = c("scenario", "row_id", "part", "pop_n", "nutrient_input")) |> 
  dplyr::group_by(scenario, part, pop_n, nutrient_input) |>
  dplyr::group_split()

#### Fit regression model ####

full_lm_df <- purrr::map_dfr(results_combined_list, function(df_temp) {
  
  df_temp_stand <- dplyr::select(df_temp, value.cv, value.prod) |> 
    dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) log(x))) # |> 
    # dplyr::mutate(dplyr::across(where(is.numeric), .fns = function(x) (x - mean(x)) / sd(x)))
  
  lm_temp <- lm(value.prod ~ value.cv, data = df_temp_stand)
  
  broom::tidy(lm_temp) |> 
    dplyr::mutate(scenario = unique(df_temp$scenario), part = unique(df_temp$part), 
                  pop_n = unique(df_temp$pop_n), nutrient_input = unique(df_temp$nutrient_input), .before = term) |> 
    dplyr::mutate(p.value.class = dplyr::case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", 
                                                   p.value < 0.05 ~ "*", p.value >= 0.05 ~ ""), 
                  r2 = summary(lm_temp)$adj.r.squared, 
                  direction = dplyr::case_when(estimate < 0.0 ~ "decrease", 
                                               estimate > 0.0 ~ "increase"))
})

#### Setup ggplot ####

size_base <- 10.0

#### Create ggplot model parameters ####

pp_combined_df <- dplyr::bind_rows(results_combined_list)

gg_pe_pp <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
  df_temp <- dplyr::filter(pp_combined_df, scenario == scenario_i, pop_n != 0) |> 
    dplyr::mutate(value.cv = log(value.cv), value.prod = log(value.prod)) # |> 
  # dplyr::group_by(nutrient_input, part) |> 
  # dplyr::mutate(value.cv = (value.cv - mean(value.cv)) / sd(value.cv), 
  #               value.prod = (value.prod - mean(value.prod)) / sd(value.prod))
  
  # ggplot
  gg_scenario <- ggplot(data = df_temp, aes(x = value.cv, y = value.prod, color = pop_n)) +
    
    # adding geoms
    geom_point(alpha = 0.1, size = 0.5) +
    geom_smooth(method = "lm", formula = "y ~ x", se = FALSE, linewidth = 0.5) +
    
    # facet grid
    facet_wrap(. ~ nutrient_input * part, nrow = 3, ncol = 3, scales = "free") +
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
    annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
    annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
    annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
    
    # set scales and labs
    scale_color_viridis_d(name = "Population size", option = "C") +
    
    # setup theme
    guides(color = guide_legend(nrow = 1)) +
    labs(x = "log(Portfolio effect)", y = expression(paste("log(Primary production [", gDW~d^-1~m^-2, "])"))) +
    theme_classic(base_size = size_base) +
    theme(strip.background = element_blank(), strip.text = element_blank(),
          axis.line = element_blank(), legend.position = "bottom")
  
  cowplot::ggdraw(gg_scenario, xlim = c(0, 1.025), ylim = c(0.0,  1.025)) +
    cowplot::draw_label("Aboveground", x = 0.2, y = 1.0, angle = 0, size = size_base * 0.65) + 
    cowplot::draw_label("Belowground", x = 0.525, y = 1.0, angle = 0, size = size_base * 0.65) + 
    cowplot::draw_label("Total", x = 0.875, y = 1.0, angle = 0, size = size_base * 0.65) + 
    cowplot::draw_label("Nutr. input: low", x = 1.0, y = 0.8, angle = 270, size = size_base * 0.65, hjust = 1) +
    cowplot::draw_label("Nutr. input: medium", x = 1.0, y = 0.575, angle = 270, size = size_base * 0.65) +
    cowplot::draw_label("Nutr. input: high", x = 1.0, y = 0.225, angle = 270, size = size_base * 0.65, hjust = 1)
  
})

#### Save ggplot ####

suppoRt::save_ggplot(plot = gg_pe_pp$noise, 
                     filename = ifelse(near, yes = "Figure-3-near.pdf", no = "Figure-3.pdf"),
                     path = "04_Figures/", width = width, height = height * 0.5,
                     units = units, dpi = dpi, overwrite = FALSE)

# #### Setup ggplot ####
# 
# color_part <- c("ag_production" = "#b6cb9d", "bg_production" = "#a6970f", "ttl_production" = "#1e4e1a")
# 
# width_pos <- 0.35
# 
# #### Create ggplot model parameters ####
# 
# gg_coef_scenario <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
#   
#   df_temp <- dplyr::filter(full_lm_df, scenario == scenario_i, term == "value.cv") 
#   
#   # ggplot
#   ggplot(data = df_temp, aes(x = pop_n, y = estimate, color = part, fill = part, group = part)) +
#     
#     # adding geoms
#     geom_hline(yintercept = 0.0, linetype = 2, color = "grey") +
#     geom_linerange(aes(ymin = 0.0, ymax = estimate), alpha = 0.75, linewidth = 0.75,
#                    position = position_dodge(width = width_pos)) +
#     geom_point(aes(shape = direction), size = 5.0,
#                position = position_dodge(width = width_pos)) +
#     geom_text(aes(label = p.value.class), vjust = 0.75, size = 2.5, color = "black",
#               position = position_dodge(width = width_pos)) +
#       
#     # facet grid
#     facet_grid(rows = dplyr::vars(nutrient_input), scales = "free",
#                labeller = labeller(nutrient_input = function(x) paste0("Nutr. input: ", x))) +
#     annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, linewidth = 1) +
#       
#     # set scales and labs
#     scale_color_manual(name = "", values = color_part) +
#     scale_fill_manual(name = "", values = color_part, 
#                       labels = c("ag_production" = "Aboveground", "bg_production" = "Belowground", 
#                                  "ttl_production" = "Total")) +
#     scale_x_discrete(limits = rev(levels(full_lm_df$pop_n))) +
#     scale_y_continuous(breaks = seq(-quantile(abs(df_temp$estimate), 0.9), quantile(abs(df_temp$estimate), 0.9), length.out = 5), 
#                        limits = function(x) c(-max(abs(df_temp$estimate)), max(abs(df_temp$estimate))), 
#                        labels = function(x) round(x, 2)) +
#     
#     # coords
#     coord_flip() +
#       
#     # setup theme
#     labs(x = "Population size", y = "Parameter estimate") +
#     guides(shape = "none", color = "none", fill = guide_legend(override.aes = list(shape = 24, color = "white"))) +
#     theme_classic(base_size = 12.0) +
#     theme(strip.background = element_blank(), strip.text = element_text(hjust = 0.5),
#           legend.position = "bottom")
#   
# })
# 
# #### Save ggplot ####
# 
# suppoRt::save_ggplot(plot = gg_coef_scenario$noise, filename = "Figure-3.pdf",
#                      path = "04_Figures/", width = width, height = height * 0.65,
#                      units = units, dpi = dpi, overwrite = FALSE)
