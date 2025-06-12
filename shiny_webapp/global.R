# global.R
# Load all necessary packages for the app

# For the Shiny framework
library(shiny)
library(shinythemes) # To make the app look nice
# library(shinyFiles) # To make file path locating easier with local runs

# For data manipulation and plotting
library(tidyverse)
library(patchwork)

# For multi-threading
library(furrr) 

# For interactive tables
library(DT)

# For bioinformatics tasks
library(Biostrings)
library(metablastr)
library(valr)

# For C++ integration
library(Rcpp)

# Custom functions for the probe design
source("../src/HCR_design_shiny.R")

# Source the C++ functions
sourceCpp("../src/thermo_calc.cpp")
sourceCpp("../src/find_optimal_probe_set_dp_cpp.cpp")