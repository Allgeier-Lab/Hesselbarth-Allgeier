##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# Purpose: Sample parameter space for simulation experiment

#### Load setup ####

source("01_Functions/setup.R")

#### Setup experiment ####

pop_n <- c(0, 8, 16, 32, 64, 128)

set.seed(42)

matrix_lhs <- lhs::improvedLHS(n = iterations, k = 2, dup = 2)

matrix_lhs[, 1] <- qunif(matrix_lhs[, 1], 0.0, 1.0) 

matrix_lhs[, 2] <- qunif(matrix_lhs[, 2], 0.0, 1.0) 

table(cut(matrix_lhs[, 1], breaks = seq(0.0, 1, 0.2)),
      cut(matrix_lhs[, 2], breaks = seq(0.0, 1, 0.2)))

experiment_df <- tibble::tibble(biotic = matrix_lhs[, 1], 
                                abiotic = matrix_lhs[, 2]) |> 
  dplyr::slice(rep(x = 1:dplyr::n(), times = length(pop_n))) |> 
  dplyr::mutate(pop_n = rep(x = pop_n, each = iterations))

nrow(experiment_df)

# experiment_df[experiment_df$pop_n == 0, "biotic"] <- 0.0

#### Check parameter space ####

dplyr::filter(experiment_df, pop_n == 8) |> 
  ggplot() + 
  geom_point(aes(x = biotic, y = abiotic)) + 
  geom_hline(yintercept = 0.0, color = "grey", linetype = 2) + geom_hline(yintercept = 1.0, color = "grey", linetype = 2) + 
  geom_vline(xintercept = 0.0, color = "grey", linetype = 2) + geom_vline(xintercept = 1.0, color = "grey", linetype = 2) + 
  labs(x = "Diversity consumer behavior", y = "Variability nutrient subsidies") +
  coord_fixed(ratio = 1) + 
  theme_classic()

#### Save experiment ####

suppoRt::save_rds(object = experiment_df, filename = "experiment-parameters.rds", 
                  path = "02_Data/", overwrite = FALSE)
