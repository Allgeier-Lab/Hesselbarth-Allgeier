##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Sample parameter space for simulation experiment

#### Load setup ####

source("05_Various/setup.R")

#### Setup experiment ####

set.seed(42)

reps <- 250

matrix_lhs <- lhs::improvedLHS(n = reps, k = 2, dup = 2)

matrix_lhs[, 1] <- qunif(matrix_lhs[, 1], 0.1, 1.0) 

matrix_lhs[, 2] <- qunif(matrix_lhs[, 2], 0.1, 1.0) 

table(cut(matrix_lhs[, 1], breaks = seq(0.1, 1, 0.1)),
      cut(matrix_lhs[, 2], breaks = seq(0.1, 1, 0.1)))

df_experiment <- tibble::tibble(biotic = matrix_lhs[, 1], 
                                abiotic = matrix_lhs[, 2]) %>% 
  dplyr::slice(rep(x = 1:dplyr::n(), times = 4)) %>% 
  dplyr::mutate(pop_n = rep(x = c(8, 16, 32, 64), each = reps))

#### Check parameter space ####

ggplot(data = df_experiment) + 
  geom_point(aes(x = biotic, y = abiotic)) + 
  geom_hline(yintercept = 0.1, color = "grey", linetype = 2) + geom_hline(yintercept = 1.0, color = "grey", linetype = 2) + 
  geom_vline(xintercept = 0.1, color = "grey", linetype = 2) + geom_vline(xintercept = 1.0, color = "grey", linetype = 2) + 
  coord_equal() + 
  theme_classic()

purrr::walk(unique(df_experiment$pop_n), function(i) {
  
  df_temp <- dplyr::filter(df_experiment, pop_n == i)
  
  tab_temp <- table(
    cut(df_temp$biotic, breaks = seq(0.1, 1, 0.1)),
    cut(df_temp$abiotic, breaks = seq(0.1, 1, 0.1))
  )
  
  if (sum(tab_temp == 0) > 0) warning(paste0(i, "; No parameter for some combination"))
  
  message(paste0("Pop = ", i, "; mean = ", round(mean(tab_temp), 2), "; sd = ", round(sd(tab_temp), 2)), 
          "; min = ", min(tab_temp), "; max = ", max(tab_temp))
  
})

#### Save experiment ####

suppoRt::save_rds(object = df_experiment, filename = "00_df_experiment.rds", 
                  path = "02_Data/", overwrite = FALSE)
