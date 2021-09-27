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

library(raster)
library(suppoRt) # remotes::install_github("mhesselbarth/suppoRt")
library(tidyverse)

library(future)
library(future.batchtools)
library(future.apply)

# Set some plotting defaults

# DINA4
units <- "mm"

width <- 210

height <- 297

# set pixels
dpi <- 900
