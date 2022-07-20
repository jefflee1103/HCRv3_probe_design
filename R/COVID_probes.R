##
# COVID probes
##

# Environment ---------------------------------------------------------------------------------

library(tidyverse)
source("./src/HCR_design.R")


# Nucleocapsid --------------------------------------------------------------------------------

## Get pre-designed sequences from SS lab, remove pre-appended B1 hairpin and attach B# of my choice

first <- read_delim(pipe("pbpaste"), col_names = "first", delim = "\t") 
second <- read_delim(pipe("pbpaste"), col_names = "second", delim = "\t")
SARS_N <- bind_cols(first, second)
SARS_N 


SARS_N <- SARS_N %>%
  mutate(first = str_sub(first, start = 1, end = 25),
         second = str_sub(second, start = 21, end = 45))
SARS_N

tmp <- attach_hcr_initiatior_standalone(SARS_N,
                                 "first",
                                 "second",
                                 "B3",
                                 "SARS2-N") %>%
  dplyr::select(5,6,7,8)


