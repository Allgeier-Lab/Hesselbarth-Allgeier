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

n <- 5

df_results <- import_data(path =  paste0("02_Data/05-variability-phase-", n, ".rds"))

#### Model dredge ####

df_dredge <- dplyr::filter(df_results, part %in% c("ag_production", "bg_production", "ttl_production"), 
                           measure %in% c("alpha", "gamma")) %>% 
  dplyr::group_by(part, measure) %>% 
  dplyr::group_split() %>% 
  purrr::map_dfr(function(df_temp) {
    
    df_temp$value.cv <- log(df_temp$value.cv) - mean(log(df_temp$value.cv))
    
    df_temp$biotic <- log(df_temp$biotic)
    df_temp$abiotic <- log(df_temp$abiotic)
    df_temp$pop_n <- as.numeric(paste(df_temp$pop_n))
    
    lm(value.cv ~ biotic + abiotic + pop_n + biotic:abiotic + biotic:pop_n + abiotic:pop_n,
       data = df_temp, na.action = "na.fail") %>%
      MuMIn::dredge() %>%
      as.data.frame() %>% 
      head(n = 1) %>% 
      dplyr::mutate(part = unique(df_temp$part), measure = unique(df_temp$measure), 
                    .before = "(Intercept)")}) %>% 
  dplyr::mutate(dplyr::across(c(3:9, 11:14), round, 3)) %>%
  dplyr::select(-c(df, logLik, AICc, weight)) %>% 
  dplyr::arrange(measure)
