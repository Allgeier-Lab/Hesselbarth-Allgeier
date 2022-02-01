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
library(terra)
library(rslurm)
library(suppoRt) # remotes::install_github("mhesselbarth/suppoRt")
library(tidyverse)
library(viridis)

#### Load data ####

default_starting <- readRDS("02_Data/default_starting.rds")

default_parameters <- readRDS("02_Data/default_parameters.rds")

#### Basic parameters ####

# set min_per_i
min_per_i <- 120

# run the model for n years
years <- 50
max_i <- (60 * 24 * 365 * years) / min_per_i

# run seagrass only 1 day
days_seagrass <- 1
seagrass_each <- (24 / (min_per_i / 60)) * days_seagrass

# save results only every m days
days_save <- 125 # which(max_i %% ((24 / (min_per_i / 60)) * (1:365)) == 0)
save_each <-  (24 / (min_per_i / 60)) * days_save

# max_i %% save_each

# set frequency of input peaks
freq_mn <- years * 1/4

# number of local metaecosystems
n <- 9

# setup extent and grain
dimensions <- c(50, 50)

grain <- 1

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

# path for rslurm
# run file.path(R.home("bin"), "Rscript") on HPC
rscript_path <- "/sw/arcts/centos7/stacks/gcc/8.2.0/R/4.1.0/lib64/R/bin/Rscript"

# check file in 99_various/hpc for raw template
sh_template <- "05_Various/rslurm_slurm.tmpl"

overwrite <- FALSE

message("\nUsing R v", stringr::str_split(rscript_path, pattern = "/", simplify = TRUE)[, 9], " on HPC")

# exclude slow nodes
exclude_nodes <- c("gl[3324-3327]", "gl[3368-3371]")

#### Remove some basic parameters ####

rm(days_save, days_seagrass)
