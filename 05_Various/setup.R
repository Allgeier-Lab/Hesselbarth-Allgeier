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

library(cowplot)
library(ggpubr)
library(magrittr)
library(Rcpp)
library(rslurm)
library(suppoRt) # remotes::install_github("mhesselbarth/suppoRt")
library(terra)
library(tidyverse)
library(vegan)

#### Load data ####

starting_list <- readRDS("02_Data/starting_list.rds")

parameters_list <- readRDS("02_Data/parameters_list.rds")

#### Basic parameters ####

# number of iterations
iterations <- 50

# set min_per_i
min_per_i <- 120

# ru<- the model for n years
years <- 50
max_i <- (60 * 24 * 365 * years) / min_per_i

# run seagrass only 1 day
days_seagrass <- 1
seagrass_each <- (24 / (min_per_i / 60)) * days_seagrass

# save results only every m days
days_save <- 25 # which(max_i %% ((24 / (min_per_i / 60)) * (1:365)) == 0)
save_each <- (24 / (min_per_i / 60)) * days_save

# max_i %% save_each

# years used to filter
years_filter <- 40

# set frequency of input peaks
freq_mn <- years * 1/4

# number of local metaecosystems
n <- 5

# setup extent and grain
dimensions <- c(50, 50)

grain <- 1

# use log size distribution (or not)
use_log <- FALSE

nutrient_input <- 4.584905e-05 # see 00_input-nutr-fish.R

#### Setup reef cells ####

reef_matrix <- matrix(data = c(-1, 0, 0, 1, 1, 0, 0, -1, 0, 0), 
                      ncol = 2, byrow = TRUE)

#### Setup treatment levels ####

amplitude_levels <- c(0.05, 0.5, 0.95)

enrichment_levels <- c(0.75, 1.0, 1.25)

#### Plotting defaults ####

# DINA4
units <- "mm"

width <- 210

height <- 297

# set pixels
dpi <- 900

# set base_size
base_size <- 12

#### HPC settings ####

account <- "jeallg0"

# path for rslurm
# run file.path(R.home("bin"), "Rscript") on HPC
rscript_path <- "/sw/arcts/centos7/stacks/gcc/8.2.0/R/4.1.0/lib64/R/bin/Rscript"

# check file in 99_various/hpc for raw template
sh_template <- "05_Various/rslurm_slurm.tmpl"

overwrite <- FALSE

message("\nUsing R v", stringr::str_split(rscript_path, pattern = "/", simplify = TRUE)[, 9], " on HPC")

# # exclude slow nodes
exclude_nodes <- "gl[3368-3371]"

#### Remove some basic parameters ####

# rm(days_save, days_seagrass)
