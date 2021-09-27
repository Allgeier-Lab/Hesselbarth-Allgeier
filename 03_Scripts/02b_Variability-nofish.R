##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

source("05_Various/setup.R")

#### Load data ####

variability_exp <- readRDS("02_Data/variability_experiment.rds")

#### Setup future plan ####

# login node -> cluster nodes -> core
login <- tweak(remote, workers = "greatlakes.arc-ts.umich.edu", user = "mhessel")

sbatch <- tweak(batchtools_slurm, template = "future_slurm.tmpl",
                resources = list(job_name = "calc-cv",
                                 log_file = "calc-cv.log",
                                 walltime = "01:00:00", # walltime <hh:mm:ss>
                                 mem_cpu  = "14G")) # memory per core in mb

plan(list(
  login,
  sbatch,
  sequential
))

### Sample CV on HPC ####

variability_nofish %<-% future.apply::future_lapply(1:nrow(variability_exp), FUN = function(i) {
  
  result %<-% {
    
    # create vector with parts to calc variability
    part <- c("bg_biomass", "ag_biomass", "bg_production", "ag_production")
    
    # create vector for lag
    lag <- c(FALSE, FALSE, TRUE, TRUE)
    
    # create filename
    file_name <- paste0("/home/mhessel/results/meta-rn_nofish_", i, ".rds")
    
    # read data
    meta_temp <- readRDS(file = file_name) 
    
    # sample variability for input
    cv_temp_in <- cbind(what = "input", part = NA,
                        meta.arrR::sample_variability(x = meta_temp$nutr_input, 
                                                      verbose = FALSE))
    
    # sample variability for output and biomass and production
    cv_temp_out <- purrr::map_dfr(seq_along(part), function(j) {
      cbind(what = "output", part = part[j],
            meta.arrR::sample_variability(x = meta_temp, what = part[j], lag = lag[j], 
                                          verbose = FALSE))
      
    })
    
    # combine to one df
    cv_combined <- rbind(cv_temp_in, cv_temp_out)
    
    # create filename
    file_name <- paste0("/home/mhessel/results/variability_nofish_", i, ".rds")
    
    # save result explicit in folder
    saveRDS(object = cv_combined, file = file_name)
    
    # only return string
    file_name
    
  }
})

#### Save data ####

# get id of all results
result_id <- list.files(path = "~/Downloads/results",
                        full.names = TRUE, pattern = "^variability_nofish_*") %>%
  stringr::str_sort(numeric = TRUE) %>% 
  map_int(function(i) stringr::str_extract(i, pattern = "[0-9]+") %>% as.integer)

# check if all ids are there
(missing_id <- which(!1:nrow(variability_exp) %in% result_id))

# read data
variability_nofish <- list.files(path = "~/Downloads/results",
                                full.names = TRUE, pattern = "^variability_nofish_*") %>%
  stringr::str_sort(numeric = TRUE) %>% 
  purrr::map_dfr(readRDS, .id = "input") %>% 
  dplyr::mutate(input = as.integer(input)) %>% 
  dplyr::left_join(y = cbind(input = 1:nrow(variability_exp), variability_exp), 
                   by = "input")

suppoRt::save_rds(object = variability_nofish, filename = "variability_nofish.rds", 
                  path = "02_Data/", overwrite = FALSE)
