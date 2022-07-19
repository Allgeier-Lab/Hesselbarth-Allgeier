##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

#### Load setup ####

source("05_Various/setup.R")

#### Adapt parameters ####

list_starting$pop_n <- 0

list_parameters$nutrients_diffusion <- 0.0

list_parameters$detritus_diffusion <- 0.0

list_parameters$detritus_fish_diffusion <- 0.0

#### Stable values ####

list_stable <- arrR::get_stable_values(bg_biomass = list_starting$bg_biomass,
                                         ag_biomass = list_starting$ag_biomass,
                                         parameters = list_parameters)

list_starting$nutrients_pool <- list_stable$nutrients_pool

list_starting$detritus_pool <- list_stable$detritus_pool

input_mn <- list_stable$nutr_input

#### Create function ####

# create globals
globals <- list(n = n, max_i = max_i, input_mn = input_mn,
                list_starting = list_starting, list_parameters = list_parameters, 
                dimensions = dimensions, grain = grain,
                amplitude = amplitude, phase = phase,
                min_per_i = min_per_i, seagrass_each = seagrass_each, save_each = save_each) 

foo <- function(freq_mn, globals) {
  
  # setup metaecosystems
  metasyst_temp <- meta.arrR::setup_meta(n = globals$n, max_i = globals$max_i, 
                                         starting_values = globals$list_starting, 
                                         parameters = globals$list_parameters, 
                                         dimensions = globals$dimensions, grain = globals$grain, 
                                         reef = NULL, verbose = FALSE)
  
  # simulate input 
  input_temp <- meta.arrR::sim_nutr_input(n = globals$n, max_i = globals$max_i,
                                          variability = c(globals$amplitude, globals$phase),
                                          input_mn = globals$input_mn, freq_mn = freq_mn, 
                                          verbose = FALSE)
  
  # run model
  result_temp <- meta.arrR::run_meta(metasyst = metasyst_temp, nutr_input = input_temp,
                                     parameters = globals$list_parameters,
                                     max_i = globals$max_i, min_per_i = globals$min_per_i, 
                                     seagrass_each = globals$seagrass_each,
                                     save_each = globals$save_each, verbose = FALSE)
  
  return(result_temp)
}

#### Submit to HPC #### 

nofish_sbatch <- rslurm::slurm_apply(f = foo, params = freq_mn_df, 
                                     globals = globals, jobname = "nofish_freq",
                                     nodes = nrow(freq_mn_df), cpus_per_node = 1, 
                                     slurm_options = list("account" = account, 
                                                          "partition" = "standard",
                                                          "time" = "00:30:00", ## hh:mm::ss
                                                          "mem-per-cpu" = "7G", 
                                                          "exclude" = "gl3324,gl3325,gl3326"),
                                     pkgs = "meta.arrR",
                                     rscript_path = rscript_path, sh_template = sh_template, 
                                     submit = FALSE)

#### Collect results ####

# nofish_freq <- rslurm::get_slurm_out(nofish_sbatch, outtype = "raw")

# rslurm::cleanup_files(nofish_sbatch)

##### Look at beta to freq ####

file_names <- list.files("_rslurm_nofish_freq/", pattern = "result_*", full.names = TRUE) %>% 
  stringr::str_sort(numeric = TRUE) 

freq_res <- purrr::map_dfr(seq_along(file_names), function(i) {
  
  message("\r> Progress: ", i, "/", length(file_names), appendLF = FALSE)
  
  x <- readRDS(file_names[i]) %>% 
    magrittr::extract2(1)
    
  beta <- meta.arrR::calc_variability(x) %>% 
    magrittr::extract2("beta")
  
  data.frame(freq = x$nutr_input$freq_mn, beta = beta)
})

ggplot(data = freq_res, aes(x = freq, y = beta)) + 
  geom_point() + geom_smooth() + 
  geom_hline(yintercept = 1.25, linetype = 2, col = "grey") + 
  geom_vline(xintercept = years * 1/4, linetype = 2, col = "grey") +
  geom_vline(xintercept = years * 1/3, linetype = 2, col = "grey") + 
  theme_classic()

j <- 25

run_temp <- list.files("_rslurm_nofish_freq/", pattern = "result_*", full.names = TRUE) %>% 
  stringr::str_sort(numeric = TRUE) %>% 
  magrittr::extract(j) %>% 
  readRDS() %>% 
  magrittr::extract2(1)

run_temp$nutr_input$freq_mn / years
plot(run_temp, summarize = TRUE)
sample_variability(run_temp, what = "ag_production", lag = TRUE)
