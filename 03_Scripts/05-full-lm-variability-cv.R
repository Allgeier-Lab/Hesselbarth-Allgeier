##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Variation partitioning noise variability

#### Load setup ####

source("05_Various/setup.R")
source("01_Functions/import_data.R")

#### Load/wrangle simulated data ####

n <- 5

results_phase_df <- import_data(path = paste0("02_Data/result-phase-", n, ".rds"))

results_noise_df <- import_data(path = paste0("02_Data/result-noise-", n, ".rds"))

results_combined_df <- dplyr::bind_rows(phase = results_phase_df, noise = results_noise_df, 
                                        .id = "scenario") %>% 
  dplyr::mutate(scenario = factor(scenario, levels = c("phase", "noise")))

#### Wrangle data ####

results_combined_list <- dplyr::filter(results_combined_df, 
                                       part %in% c("ag_production", "bg_production", "ttl_production"),
                                       measure %in% c("alpha", "gamma")) %>% 
  dplyr::group_by(scenario, part, measure) %>% 
  dplyr::group_split() %>% 
  purrr::map(function(temp_df) {
   
    dplyr::select(temp_df, value.cv, value.cv.sd, value.cv.mn, biotic, abiotic, pop_n) %>% 
      dplyr::mutate(pop_n = as.numeric(paste(pop_n))) %>%
      dplyr::mutate(across(.fns = function(x) log(x))) %>%
      dplyr::mutate(across(.fns = function(x) (x - mean(x)) / sd(x))) %>% 
      dplyr::mutate(scenario = unique(temp_df$scenario), part = unique(temp_df$part), measure = unique(temp_df$measure), 
                    .before = "value.cv") 
    
})

names(results_combined_list) <- purrr::map_chr(results_combined_list, function(i) {
  paste(unique(i$scenario), unique(i$part), unique(i$measure), sep = "-")
})

#### Run linear regression model ####

estimators_df <- purrr::map_dfr(results_combined_list, function(temp_df) {
    
    lm(value.cv ~ biotic + abiotic + pop_n + biotic:abiotic + biotic:pop_n + abiotic:pop_n,
       data = temp_df, na.action = "na.fail") %>% 
    broom::tidy()}, .id = "id") %>% 
  tidyr::separate(col = id, into = c("scenario", "part", "measure"), sep = "-") %>% 
  dplyr::mutate(estimate = dplyr::case_when(p.value < 0.05 ~ estimate,
                                            p.value > 0.05 ~ as.numeric(NA))) %>% 
  dplyr::select(-c(std.error, statistic, p.value)) %>% 
  tidyr::pivot_wider(names_from = term, values_from = estimate) %>% 
  dplyr::arrange(scenario, part, measure)

# clipr::write_clip(estimators_df)

#### Model dredge ####

dredge_df <- purrr::map_dfr(results_combined_list, function(temp_df) {
  
  lm(value.cv ~ biotic + abiotic + pop_n + biotic:abiotic + biotic:pop_n + abiotic:pop_n,
     data = temp_df, na.action = "na.fail") %>% 
    MuMIn::dredge() %>% 
    as.data.frame() %>% 
    head(n = 1)}, .id = "id") %>% 
  tidyr::separate(col = id, into = c("scenario", "part", "measure"), sep = "-") %>%
  dplyr::arrange(scenario, part, measure) %>% 
  dplyr::select(-c(df, logLik, AICc, delta, weight)) %>% 
  tibble::as_tibble()

# clipr::write_clip(dredge_df)
