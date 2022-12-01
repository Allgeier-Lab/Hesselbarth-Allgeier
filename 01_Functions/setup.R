##--------------------------------------------##
##    Author: Maximilian H.K. Hesselbarth     ##
##    Coastal Ecology and Conservation Lab    ##
##    University of Michigan                  ##
##    mhessel@umich.edu                       ##
##    www.github.com/mhesselbarth             ##
##--------------------------------------------##

# load packages #
library(meta.arrR) # remotes::install_github("Allgeier-Lab/meta.arrR", ref = "development")
library(arrR)

library(broom)
library(cowplot)
library(ggpmisc)
library(lhs)
library(MetBrewer)
library(MuMIn)
library(relaimpo)
library(rslurm)
library(suppoRt) # remotes::install_github("mhesselbarth/suppoRt")
library(tidyverse)

#### Load data ####

starting_values_list <- readRDS("02_Data/starting-values.rds")

parameters_list <- readRDS("02_Data/parameter-values.rds")

#### Basic parameters ####

# number of iterations
iterations <- 50

# set min_per_i
min_per_i <- 120

# the model for n years
years <- 50
max_i <- (60 * 24 * 365 * years) / min_per_i

# set frequency to once a year
frequency <- years # / 2

# run seagrass only 1 day
days_seagrass <- 1
seagrass_each <- (24 / (min_per_i / 60)) * days_seagrass

# save results only every m days
days_save <- 25 # which(max_i %% ((24 / (min_per_i / 60)) * (1:365)) == 0)
save_each <- (24 / (min_per_i / 60)) * days_save

# max_i %% save_each

# years used to filter
years_filter <- 40

# setup extent and grain
dimensions <- c(50, 50)

grain <- 1

# number of ecosystems
n <- 5

#### Setup reef cells ####

reef_matrix <- matrix(data = c(-1, 0, 0, 1, 1, 0, 0, -1, 0, 0), 
                      ncol = 2, byrow = TRUE)

#### ####

threshold_abundance <- 135

#### Plotting defaults ####

# DINA4
units <- "mm"

width <- 210

height <- 297

# set pixels
dpi <- 900

extension <- ".pdf"

#### HPC settings ####

account <- "jeallg0"

# path for rslurm
# run file.path(R.home("bin"), "Rscript") on HPC
rscript_path <- "/sw/pkgs/arc/stacks/gcc/10.3.0/R/4.2.0/lib64/R/bin/Rscript"

message("\nUsing R v", stringr::str_split(rscript_path, pattern = "/", simplify = TRUE)[, 9], " on HPC")

# # exclude slow nodes
# exclude_nodes <- "gl[3368-3371]"
# exclude_nodes <- "gl3069"

#### Remove some basic parameters ####

# rm(days_save, days_seagrass)
