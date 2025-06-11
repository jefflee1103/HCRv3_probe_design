# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Design FLARIMv2 probes
# https://www.biorxiv.org/content/10.1101/2021.08.13.456301v1.full.pdf
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Environment -------------------------------------------------------------

library(tidyverse)

# Initialtion halves that should be attached to the end of ribosome probes
initiator_I1_a <- "gAggAgggCAgCAAACgg"
initiator_I2_a <- "CCTCgTAAATCCTCATCA"
initiator_I3_a <- "gTCCCTgCCTCTATATCT"
initiator_I4_a <- "CCTCAACCTACCTCCAAC"
initiator_I5_a <- "CTCACTCCCAATCTCTAT"

# Create 18S FLARIMv2 HCRB1 probes 
richer_HCRB3 <- read_csv("./data/FLARIMv2/18S_FLARIMv2_sB3.csv")

FLARIM_18S_HCRB1 <- richer_HCRB3 %>%
  mutate(sequence = str_sub(sequence, start = 1, end = -19)) %>%
  mutate(sequence = paste0(sequence, initiator_I1_a)) %>%
  mutate(name = str_replace(name, "sB3", "sB1"))

FLARIM_18S_HCRB1

write_csv(FLARIM_18S_HCRB1, "./data/FLARIMv2/18S_FLARIMv2_sB1.csv")

# Create 18S FLARIMv2 HCRB2 probes 
richer_HCRB3 <- read_csv("./data/FLARIMv2/18S_FLARIMv2_sB3.csv")

FLARIM_18S_HCRB2 <- richer_HCRB3 %>%
  mutate(sequence = str_sub(sequence, start = 1, end = -19)) %>%
  mutate(sequence = paste0(sequence, initiator_I2_a)) %>%
  mutate(name = str_replace(name, "sB3", "sB2"))

FLARIM_18S_HCRB2

write_csv(FLARIM_18S_HCRB2, "./data/FLARIMv2/18S_FLARIMv2_sB2.csv")

# Create 18S FLARIMv2 HCRB4 probes 
richer_HCRB3 <- read_csv("./data/FLARIMv2/18S_FLARIMv2_sB3.csv")

FLARIM_18S_HCRB4 <- richer_HCRB3 %>%
  mutate(sequence = str_sub(sequence, start = 1, end = -19)) %>%
  mutate(sequence = paste0(sequence, initiator_I4_a)) %>%
  mutate(name = str_replace(name, "sB3", "sB4"))

FLARIM_18S_HCRB4

write_csv(FLARIM_18S_HCRB4, "./data/FLARIMv2/18S_FLARIMv2_sB4.csv")
