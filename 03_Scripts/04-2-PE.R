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

results_phase_df <- import_cv(path = "02_Data/result-phase.rds")

results_noise_df <- import_cv(path = "02_Data/result-noise.rds")

results_combined_df <- dplyr::bind_rows(phase = results_phase_df, noise = results_noise_df, 
                                        .id = "scenario") |> 
  dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise")), 
                type = dplyr::case_when(pop_n == 0 ~ "a", abiotic == 0.0 ~ "b", TRUE ~ "ab"), 
                type = factor(type, levels = c("a", "ab", "b"), 
                              labels = c("Abiotic subsidies only", "Abiotic and connectivity", 
                                         "Consumer connectivity only")))

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

beta_df <- dplyr::left_join(x = results_combined_df, y = mortality_combined_df, 
                                        by = c("scenario", "row_id", "pop_n", "nutrient_input")) |> 
  dplyr::mutate(include = dplyr::case_when(pop_n == 128 & nutrient_input == "low" &
                                             died_total > threshold_mort ~ "no", 
                                           TRUE ~ "yes")) |> 
  dplyr::filter(include == "yes", measure == "beta")


#### Create ggplot ####

size_base <- 10.0
size_text <- 2.0
size_point <- 2.5

width_doge <- 0.5

gg_pe <- purrr::map(c(phase = "phase", noise = "noise"), function(scenario_i) {
  
  purrr::map(c(ag = "ag_production", bg =  "bg_production", ttl = "ttl_production"), function(part_i) {
    
    beta_temp <- dplyr::filter(beta_df, scenario == scenario_i, part == part_i)
    
    label_temp <- dplyr::group_by(beta_temp, nutrient_input) |> 
      dplyr::group_split() |> 
      purrr::map_dfr(function(i) {
        anova <- aov(value.cv ~ type, data = i)
        
        tukey <- TukeyHSD(anova)
        
        label_letters <- multcompView::multcompLetters4(anova, tukey)
        
        data.frame(nutrient_input = unique(i$nutrient_input), type = names(label_letters$type$monospacedLetters),
                   letter = unname(label_letters$type$monospacedLetters))
        
      })
    
    ggplot(beta_temp, aes(x = type, y = value.cv)) + 
      
      # geoms
      geom_boxplot(outlier.shape = NA) +
      geom_hline(yintercept = 1.0, linetype = 2, color = "grey") +
      geom_text(data = label_temp, aes(x = type, y = max(beta_temp$value.cv) * 0.85, label = letter)) +
      
      # scales
      coord_cartesian(ylim = c(1.0, max(beta_temp$value.cv) * 0.85)) +
      
      # facet
      facet_grid(rows = vars(nutrient_input),
                 labeller = labeller(nutrient_input = function(x) paste("Nutr. input:", x))) +
      
      # themes
      labs(x = "", y = "Portfolio effect") +
      guides(fill = guide_legend(nrow = 1)) +
      theme_classic(base_size = 10) +
      theme(strip.background = element_blank(), strip.text = element_text(size = size_base * 0.65),
            axis.line = element_blank(), panel.border = element_rect(linewidth = 0.5, fill = NA), 
            legend.position = "bottom")
  })
})

#### save ggplot ####
suppoRt::save_ggplot(plot = gg_pe$noise$ag, filename = "Figure-A5.pdf",
                     path = "04_Figures/Appendix/", width = width, height = height * 0.5,
                     units = units, dpi = dpi, overwrite = FALSE)
