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

reps <- 100

matrix_lhs <- lhs::improvedLHS(n = reps, k = 2, dup = 2)

matrix_lhs[, 1] <- qunif(matrix_lhs[, 1], 0.1, 1.0) 

matrix_lhs[, 2] <- qunif(matrix_lhs[, 2], 0.1, 1.0) 

table(cut(matrix_lhs[, 1], breaks = seq(0.1, 1, 0.1)),
      cut(matrix_lhs[, 2], breaks = seq(0.1, 1, 0.1)))

experiment_df <- tibble::tibble(biotic = matrix_lhs[, 1], 
                                abiotic = matrix_lhs[, 2]) |> 
  dplyr::slice(rep(x = 1:dplyr::n(), times = 5)) |> 
  dplyr::mutate(pop_n = rep(x = c(8, 16, 32, 64, 128), each = reps))

nrow(experiment_df)

#### Check parameter space ####

ggplot(data = experiment_df) + 
  geom_point(aes(x = biotic, y = abiotic)) + 
  geom_hline(yintercept = 0.1, color = "grey", linetype = 2) + geom_hline(yintercept = 1.0, color = "grey", linetype = 2) + 
  geom_vline(xintercept = 0.1, color = "grey", linetype = 2) + geom_vline(xintercept = 1.0, color = "grey", linetype = 2) + 
  coord_equal() + 
  theme_classic()

purrr::walk(unique(experiment_df$pop_n), function(i) {
  
  temp_df <- dplyr::filter(experiment_df, pop_n == i)
  
  tab_temp <- table(
    cut(temp_df$biotic, breaks = seq(0.1, 1, 0.2)),
    cut(temp_df$abiotic, breaks = seq(0.1, 1, 0.2))
  )
  
  if (sum(tab_temp == 0) > 0) warning(paste0(i, "; No parameter for some combination"))
  
  message(paste0("Pop = ", i, "; mean = ", round(mean(tab_temp), 2), "; sd = ", round(sd(tab_temp), 2)), 
          "; min = ", min(tab_temp), "; max = ", max(tab_temp))
  
})

#### Save experiment ####

suppoRt::save_rds(object = experiment_df, filename = "experiment-parameters.rds", 
                  path = "02_Data/", overwrite = FALSE)
