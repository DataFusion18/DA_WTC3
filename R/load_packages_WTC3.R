# This script loads necessary packages
##############################
# Load packages
library("mvtnorm") # Creates candidate parameter vector as a multivariate normal jump away from the current candidate
library("reshape2")
library("ggplot2")
library("lubridate")
library("rio")
library("dplyr")
library("zoo")
library("doBy")
library("corrplot")
library("png")
library("grid")
library("gridExtra")
library("plyr")
library("car")
library("lattice")
library("visreg")
library("data.table")
library("grDevices")
library(devtools)
install_bitbucket("remkoduursma/plotBy")

#- load the custom functions that do most of the heavy lifting
source("R/stat_smooth_func.R")


