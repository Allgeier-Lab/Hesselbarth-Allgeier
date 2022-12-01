##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Create abundance vs mort figure

#### Load setup ####

source("01_Functions/setup.R")
source("01_Functions/import-mortality.R")
source("01_Functions/import-abundance.R")

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

#### Import abundance data ####

abundance_phase_df <- import_abundance(path = "02_Data/result-phase.rds") |> 
  dplyr::filter(pop_n > 0)

abundance_noise_df <- import_abundance(path = "02_Data/result-noise.rds") |> 
  dplyr::filter(pop_n > 0)

abundance_combined_df <- dplyr::bind_rows(phase = abundance_phase_df, noise = abundance_noise_df,
                                          .id = "scenario") |>
  dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise"))) |> 
  dplyr::group_by(scenario, row_id, pop_n, nutrient_input) |>
  dplyr::summarise(abundance_max = max(mean), .groups = "drop")

#### Combine to one data.frame ####

result_total_df <- dplyr::left_join(x = abundance_combined_df, y = mortality_combined_df, 
                                    by = c("scenario", "row_id", "pop_n", "nutrient_input"))

#### ####

# ggplot
ggplot(result_total_df, aes(x = abundance_max, y = died_total))  + 
  
  #geoms 
  geom_point(pch = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "#6ad5e8") +
  
  # facet
  # facet_grid(rows = vars(nutrient_input), cols = vars(pop_n), scales = "free") + 
  facet_wrap(. ~ nutrient_input + pop_n, scales = "free", nrow = 3, ncol = 5, 
             labeller = labeller(nutrient_input = function(x) paste0("Nutr. input: ", x), 
                                 pop_n = function(x) paste0("Pop. size: ", x))) +
  
  # themes
  labs(x = "Maximum abundance", y = "Maximum total mortality") + 
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(hjust = 0),
        axis.line = element_blank(), panel.border = element_rect(linewidth = 0.5, fill = NA))

threshold <- 135

(gg_128 <- dplyr::filter(result_total_df, pop_n == 128, nutrient_input == "low") |>
  dplyr::mutate(group = dplyr::case_when(abundance_max <= threshold ~ "low", TRUE ~ "high"),
                group = factor(group, levels = c("low", "high"))) |>
  
  # ggplot
  ggplot(aes(x = abundance_max, y = died_total))  + 
  
  # geoms 
  geom_point(pch = 1, alpha = 0.75) +
  geom_smooth(method = "lm", se = FALSE, color = "#6ad5e8") +
  # geom_vline(xintercept = threshold, linetype = 2, color = "#f7aa58") +

  # facet_wrap
  facet_wrap(. ~ group, scales = "free_x") +
  
  # themes
  labs(x = "Maximum abundance", y = "Maximum total mortality") + 
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_blank(),
        axis.line = element_blank(), panel.border = element_rect(linewidth = 0.5, fill = NA)))

#### Save ggplot ####

suppoRt::save_ggplot(plot = gg_128, filename = "Figure-A1.pdf",
                     path = "04_Figures/Appendix/", width = width, height = height * 0.35,
                     units = units, dpi = dpi, overwrite = FALSE)
