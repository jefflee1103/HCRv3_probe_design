# global.R
# Load all necessary packages for the app

# For the Shiny framework
library(shiny)
library(shinythemes) # To make the app look nice

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

# Custom functions for the probe design
source("R/HCR_design_shiny.R")

