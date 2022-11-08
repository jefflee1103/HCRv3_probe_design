##
# COVID probes
##

# Environment ---------------------------------------------------------------------------------

library(tidyverse)
source("./R/HCR_design.R")


# Nucleocapsid --------------------------------------------------------------------------------

## Get pre-designed sequences from SS lab, remove pre-appended B1 hairpin and attach B# of my choice

first <- read_delim(pipe("pbpaste"), col_names = "first", delim = "\t") 
second <- read_delim(pipe("pbpaste"), col_names = "second", delim = "\t")
probe_set <- bind_cols(first, second)
probe_set 
probe_set <- probe_set %>%
  mutate(first = str_sub(first, start = 1, end = 25),
         second = str_sub(second, start = 21, end = 45))
probe_set

output <- attach_hcr_initiatior_standalone(probe_set,
                                 "first",
                                 "second",
                                 "B1",
                                 "SARS2-N") %>%
  dplyr::select(5,6,7,8)

output

write_csv(output, "~/Desktop/SARS2-N_HCRB1_probes-to-order.csv")


# ORF1a -------------------------------------------------------------------

## Get pre-designed sequences from SS lab, remove pre-appended B1 hairpin and attach B# of my choice

first <- read_delim(pipe("pbpaste"), col_names = "first", delim = "\t") 
second <- read_delim(pipe("pbpaste"), col_names = "second", delim = "\t")
probe_set <- bind_cols(first, second)
probe_set 
probe_set <- probe_set %>%
  mutate(first = str_sub(first, start = 1, end = 25),
         second = str_sub(second, start = 21, end = 45))
probe_set

output <- attach_hcr_initiatior_standalone(probe_set,
                                           "first",
                                           "second",
                                           "B3",
                                           "SARS2-orf1a") %>%
  dplyr::select(5,6,7,8)

output

write_csv(output, "~/Desktop/SARS2-orf1a_HCRB3_probes-to-order.csv")

