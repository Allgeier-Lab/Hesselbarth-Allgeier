##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Create figure of connectivity for increasing variability

#### Load setup ####

source("05_Various/setup.R")

#### Adapt parameters ####

# number of local metaecosystems
n <- 5 # 5 9

# the model for n years
years <- 10
max_i <- (60 * 24 * 365 * years) / min_per_i

starting_values_list$pop_n <- 16

#### Stable values #### 

stable_values_list <- arrR::get_req_nutrients(bg_biomass = starting_values_list$bg_biomass,
                                              ag_biomass = starting_values_list$ag_biomass,
                                              parameters = parameters_list)

starting_values_list$nutrients_pool <- stable_values_list$nutrients_pool

starting_values_list$detritus_pool <- stable_values_list$detritus_pool

#### Setup experiment ####

experiment_df <- readRDS("02_Data/experiment-parameters.rds")

#### Run model ####

result_list <- purrr::map(c(low = min(experiment_df$biotic), medium = quantile(x = experiment_df$biotic, probs = 0.25),
                            high = max(experiment_df$biotic)), function(connect_i) {
  
  # update move meta_sd parameters
  parameters_list$move_meta_sd <- connect_i
  
  # setup metaecosystems
  metasyst_temp <- meta.arrR::setup_meta(n = n, max_i = max_i, reef = reef_matrix,
                                         starting_values = starting_values_list, parameters = parameters_list,
                                         dimensions = dimensions, grain = grain, 
                                         use_log = FALSE, verbose = FALSE)
  
  # run model
  result_temp <- meta.arrR::run_simulation_meta(metasyst = metasyst_temp, parameters = parameters_list,
                                                movement = "behav", max_i = max_i, min_per_i = min_per_i,
                                                seagrass_each = seagrass_each, save_each = save_each) # |> 
    # meta.arrR::filter_meta(filter = c((max_i / years) * years_filter, max_i), 
    #                        reset = TRUE, verbose = FALSE)
  
  abundance_df <- meta.arrR::get_abundance(result_temp)
  
  # get moved counts
  moved_df <- dplyr::bind_rows(result_temp$fishpop) |>
    dplyr::select(timestep, id, moved) |> 
    dplyr::arrange(id)
  
  return(list(abundance = abundance_df, moved = moved_df))
  })

# # combine to data.frames
# abundance_df <- dplyr::bind_rows(low = result_list$low$abundance, high = result_list$high$abundance, 
#                                  .id = "variability")

moved_df <- dplyr::bind_rows(low = result_list$low$moved, medium = result_list$medium$moved,
                             high = result_list$high$moved,.id = "variability")

#### Summarize data ####

moved_sum_df <- dplyr::group_by(moved_df, variability, timestep) |> 
  dplyr::summarise(moved_mn = mean(moved), moved_low = mean(moved) - sd(moved),
                   moved_hi = mean(moved) + sd(moved),.groups = "drop") |> 
  dplyr::mutate(variability = factor(variability, levels = c("low", "medium", "high")))

#### Setup plots ####

breaks_x <- seq(from = 0, to = max_i, by = max_i / years)
# labels_x <- paste(seq(from = 0, to = years / filter_factor, by = 1), "years")

label_size <- 7.5

line_size <- 0.25

line_color <- c(low = "#1e466e", medium = "#ffd06e", high = "#e76254")

base_size <- 7.5

extension <- ".png"

#### Create ggplot ####

gg_moved <- ggplot(data = moved_sum_df, aes(x = timestep, y = moved_mn, color = variability)) + 
  
  # adding geoms
  # geom_ribbon(aes(ymin = moved_low, ymax = moved_hi, fill = variability),
  #             color = NA, alpha = 0.1) +
  geom_line(size = line_size) +
  
  # setup scales
  scale_fill_manual(values = line_color) +
  scale_color_manual(values = line_color) +
  scale_x_continuous(breaks = breaks_x) +
  
  # themes and stuff
  labs(x = "Time", y = "Connectivity") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "none", axis.text = element_blank())

#### Save ggplot ####

suppoRt::save_ggplot(plot = gg_moved, filename = paste0("Figure-1-connect", extension),
                     path = "04_Figures/", width = 100, height = 35,
                     units = units, dpi = dpi, overwrite = T)
