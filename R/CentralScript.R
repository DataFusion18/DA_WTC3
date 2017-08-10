#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#- Central analysis script for DA with WTC3 experiment manuscript, entitled
#  "".

#  The idea is to keep this script nice and tidy, but reproducibly do all the
#  analysis and make all of the figures for the manuscript. Raw and processed data will be 
#  placed in the "raw_data" and "processed_data" folders respectively, while figures 
#  and tables will be placed in the "output" folder.
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# Clear the workspace (if needed)
rm(list=ls())
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
#- Load the packages and custom functions that do all the work.
#  This will install a bunch of required libraries, including some non-standard stuff.
source("R/load_packages_WTC3.R")
# source("R/loadLibraries_WTC3.R")

#- load the custom analysis and plotting functions that do all of the actual work
# source("R/functions_WTC3.R")
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
# #- Download data files for Sink limited pot experiment. This downloads the zipfile from figshare
# download.file("https://ndownloader.figshare.com/files/8724376", "raw_data.zip", mode="wb")
# # Extract data to different folders.
# unzip("raw_data.zip")
# 
# #- Download data files for WTC3 experiment. This downloads the zipfile from figshare
# download.file("https://ndownloader.figshare.com/files/4857112", "data.zip", mode="wb")
# # Extract data to different folders.
# unzip("data.zip")
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
#- This script imports and processes the raw Sink limited container volume experiment data 
#  to model the carbon pools and fluxes using MCMC
source("R/initial_data_processing.R")



