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

extension <- ".pdf"

#### Load/wrangle simulated data ####

amplitude <- "095"

file_path <- paste0("02_Data/05-variability-phase-", amplitude, ".rds")

df_results <- import_data(path = file_path)

#### Model dredge ####

df_dredge <- dplyr::group_by(df_results, part, measure) %>% 
  dplyr::group_split() %>% 
  purrr::map_dfr(function(df_temp) {
    
    lm(value.cv ~ as.numeric(pop_n) + move_meta_sd + phase_sd + as.numeric(pop_n):move_meta_sd + 
         as.numeric(pop_n):phase_sd + move_meta_sd:phase_sd, data = df_temp, na.action = "na.fail") %>%
      MuMIn::dredge() %>%
      as.data.frame() %>% 
      head(n = 3) %>% 
      dplyr::mutate(part = unique(df_temp$part), measure = unique(df_temp$measure), 
                    .before = "(Intercept)")
  })

# dplyr::mutate(df_dredge, dplyr::across(c(3:9, 11:14), round, 3)) %>%
#   dplyr::select(-c(df, logLik, AICc, weight)) %>% View
