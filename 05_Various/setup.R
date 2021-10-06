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
library(raster)
library(rslurm)
library(suppoRt) # remotes::install_github("mhesselbarth/suppoRt")
library(tidyverse)

# library(future)
# library(future.batchtools)
# library(future.apply)

# Set some plotting defaults

# DINA4
units <- "mm"

width <- 210

height <- 297

# set pixels
dpi <- 900

# path for rslurm
# run file.path(R.home("bin"), "Rscript") on HPC
rscript_path <- "/sw/arcts/centos7/stacks/gcc/8.2.0/R/4.0.2/lib64/R/bin/Rscript"

# check file in 99_various/hpc for raw template
sh_template <- "05_Various/rslurm_slurm.tmpl"
