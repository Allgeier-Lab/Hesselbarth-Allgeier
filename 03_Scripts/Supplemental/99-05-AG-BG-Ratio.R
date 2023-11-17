##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: 

#### Load setup ####

source("01_Functions/setup.R")
source("01_Functions/import-cv.R")
source("01_Functions/import-mortality.R")

#### Load/wrangle simulated data ####

results_noise_df <- import_cv(path = "02_Data/result-noise.rds")

#### Load mortality data ####

mortality_noise_df <- import_mortality(path = "02_Data/result-noise.rds") |> 
  dplyr::filter(pop_n != 0) |> 
  dplyr::group_by(row_id, pop_n, nutrient_input) |>
  dplyr::summarise(died_total = mean(died_total), .groups = "drop")

#### Filter abundances/mortality #### 

results_final_df <- dplyr::left_join(x = results_noise_df, y = mortality_noise_df, 
                                     by = c("row_id", "pop_n", "nutrient_input")) |> 
  dplyr::mutate(include = dplyr::case_when(died_total > threshold_mort ~ "no", TRUE ~ "yes")) |> 
  dplyr::filter(include == "yes") |>
  dplyr::mutate(treatment = dplyr::case_when(abiotic != 0.0 & biotic == 0.0 ~ "subsidies", 
                                             abiotic == 0.0 & biotic != 0.0 ~ "connectivity", 
                                             abiotic != 0.0 & biotic != 0.0 ~ "combined")) |>
  dplyr::mutate(part = factor(part, labels = c("ag_production" = "Aboveground", "bg_production" = "Belowground", 
                                               "ttl_productio" = "Total")),
                treatment = factor(treatment, levels = c("subsidies", "connectivity", "combined")))

#### Plot Total vs AG/BG ####

base_size <- 13.5

df_pp <- dplyr::filter(results_final_df, measure == "gamma", include == "yes", treatment == "combined") |> 
  dplyr::select(row_id, part, value.prod, pop_n, nutrient_input) |> 
  dplyr::filter(part != "Total") |> 
  dplyr::group_by(nutrient_input, pop_n, part) |> 
  dplyr::summarise(value.prod = mean(value.prod), .groups = "drop")
  
gg_prod <- ggplot(data = df_pp) +
  
  # adding col geom
  geom_col(aes(x = pop_n, y = value.prod, fill = part), position = position_stack()) +
  
  # facet wrap by treatments
  facet_wrap(. ~ nutrient_input, labeller = labeller(nutrient_input = function(x) paste("Nutr. enrich.:", x))) +
  
  # adding box
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf) +
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
  
  # set scales
  scale_fill_manual(name = "Components", values = c(Aboveground = "#44AA99", Belowground = "#117733")) +
  
  # set themes
  labs(y = expression(paste("Total primary production [", gDW~d^-1~m^-2, "]")), 
       x = "Population size") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom", axis.line = element_blank(),
        strip.background = element_blank(), strip.text = element_text(hjust = 0))

#### Production values ####

dplyr::filter(results_final_df, measure == "alpha", include == "yes", treatment == "combined") |> 
  dplyr::select(row_id, part, value.prod, pop_n, nutrient_input) |> 
  dplyr::filter(part != "Total") |> 
  dplyr::group_by(nutrient_input, part) |>
  dplyr::summarise(prod.mn = mean(value.prod), prod.sd = sd(value.prod), n= dplyr::n(),
                   .groups = "drop")

#### Save figures #### 

suppoRt::save_ggplot(plot = gg_prod, filename = "Figure-S5.png",
                     path = "04_Figures/Supplemental/", width = width, height = height * 0.35,
                     units = units, dpi = dpi, overwrite = FALSE)

suppoRt::save_ggplot(plot = gg_prod, filename = "Figure-S5.pdf",
                     path = "04_Figures/Supplemental/", width = width, height = height * 0.35,
                     units = units, dpi = dpi, overwrite = FALSE)
