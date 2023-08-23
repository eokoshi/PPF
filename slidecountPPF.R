## Header --------------------------- 
##
## Script name: slidecountPPF.R
##
## Purpose of script: Count the number of slides per case and calculate the average; PPF paper
##
## Author: Ethan N. Okoshi
##
## Date Created: 2023-03-20
##
## Copyright (c) Ethan Okoshi, 2023
## Email: ethanokoshi@gmail.com
##
## Notes:
##   
## set working directory

# setwd()

## set options

options(scipen = 6, digits = 4) # I prefer to view outputs in non-scientific notation
Sys.setenv(LANG = "EN")
extrafont::loadfonts(device = "win")

## load functions from another file

# source() 

## load data

input_data_path <- "input/PPF for USCAP230309_R.xlsx"
library(readxl)
slides <- read_excel(input_data_path,sheet = "Sheet1")

## load packages -------------

require(tidyverse)

## End Header ---------

count_slides <- function(sheet) {
  mean_slides <- sheet %>% 
    select(contains("true")) %>%
    mutate(slidecount = rowSums(!is.na(.)), .keep = "none",
           S5 = !is.na(.$true_UIP_rate_S5_ishi_),
           S8 = !is.na(.$true_UIP_rate_S8_ishi),
           S9 = !is.na(.$true_UIP_rate_S9_ishi),) %>%
    mutate(across(S5:S9, as.numeric)) %>%
    bind_cols(sheet[1:2],.)
  return(mean_slides)
}
slidecount_tibble <- count_slides(slides)
mean(slidecount_tibble$slidecount)

count(slidecount_tibble, slidecount)

filter(slidecount_tibble, slidecount == 2) %>% #whats missing? two slides
  mutate(leftout = case_when(S5 == 0 ~ "S5",
                             S8 == 0 ~ "S8",
                             S9 == 0 ~ "S9")) %>%
  count(leftout)

filter(slidecount_tibble, slidecount == 1) %>% #whats there? single slide
  mutate(singleslide = case_when(S5 == 1 ~ "S5",
                             S8 == 1 ~ "S8",
                             S9 == 1 ~ "S9")) %>%
  count(singleslide)

chisq.test(sheet$focalUIP_consensus, slidecount_tibble$slidecount)$expected
fisher.test(sheet$focalUIP_consensus, slidecount_tibble$slidecount)
