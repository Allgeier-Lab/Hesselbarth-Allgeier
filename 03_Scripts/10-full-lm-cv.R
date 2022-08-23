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
source("05_Various/import_data.R")

#### Load/wrangle simulated data ####

n <- 5

phase_df <- import_data(path = paste0("02_Data/05-variability-phase-", n, ".rds"))

noise_df <- import_data(path =  paste0("02_Data/05-variability-noise-", n, ".rds"))

combined_df <- dplyr::bind_rows(phase = phase_df, noise = noise_df, .id = "scenario")

combined_list <- dplyr::filter(combined_df, part %in% c("ag_production", "bg_production", "ttl_production"),
                               measure %in% c("alpha", "gamma")) %>% 
  dplyr::group_by(scenario, part, measure) %>% 
  dplyr::group_split() %>% 
  purrr::map(function(temp_df) {
   
    temp_df$value.cv <- log(temp_df$value.cv) - mean(log(temp_df$value.cv))
    
    temp_df$biotic <- log(temp_df$biotic)
    temp_df$abiotic <- log(temp_df$abiotic)
    temp_df$pop_n <- as.numeric(paste(temp_df$pop_n))
    
    return(temp_df)
    
})

names(combined_list) <- purrr::map_chr(combined_list, function(i) {
  paste(unique(i$scenario), unique(i$part), unique(i$measure), sep = "-")
})

#### Run linear regression model ####

estimators_df <- purrr::map_dfr(combined_list, function(temp_df) {
    
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

dredge_df <- purrr::map_dfr(combined_list, function(temp_df) {
  
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
