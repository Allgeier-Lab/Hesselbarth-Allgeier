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

df_phase <- import_data(path = paste0("02_Data/05-variability-phase-", amplitude, ".rds"))

df_noise <- import_data(path =  paste0("02_Data/05-variability-noise-", amplitude, ".rds"))

df_total <- dplyr::bind_rows(phase = df_phase, noise = df_noise, .id = "scenario")

#### Fit regression model ####

df_regression <- dplyr::filter(df_total, part %in% c("ag_production", "bg_production", "ttl_production"), 
                               measure %in% c("alpha", "gamma")) %>%
  dplyr::group_by(scenario, part, measure) %>%
  dplyr::group_split() %>%
  purrr::map_dfr(function(df_temp) {
    
    df_temp$value.prod <- log(df_temp$value.prod) - mean(log(df_temp$value.prod))
    
    df_temp$biotic <- log(df_temp$biotic)
    df_temp$abiotic <- log(df_temp$abiotic)
    df_temp$pop_n <- as.numeric(paste(df_temp$pop_n))
    
    lm_temp <- lm(value.prod ~ biotic + abiotic + pop_n + 
         biotic:abiotic + biotic:pop_n + abiotic:pop_n, data = df_temp)
    
    broom::tidy(lm_temp) %>%
      dplyr::mutate(scenario = unique(df_temp$scenario), part = unique(df_temp$part), 
                    measure = unique(df_temp$measure), .before = term) %>% 
      dplyr::mutate(term = dplyr::case_when(term == "(Intercept)" ~ "intercept", TRUE ~ term), 
                    p.value = dplyr::case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**",
                                               p.value < 0.05 ~  "*", p.value >= 0.05 ~ ""),
                    direction = dplyr::case_when(term != "intercept" & estimate < 0.0 ~ "negative", 
                                                 term != "intercept" & estimate > 0.0 ~ "positive", 
                                                 term == "intercept" ~ as.character(NA)))
  })

