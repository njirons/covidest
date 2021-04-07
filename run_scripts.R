# setwd("~/Documents/research/raftery/final_results")

#setwd("~/final_results"))  ## madrid servers
setwd("~/20210223results")  ## madrid servers
#load in covidtracking data
library(stats4)
library(readr)
library(tidyverse)
library(coda)
library(lubridate)
library(Epi)
library(latex2exp)
library(rstan)
#options(mc.cores = parallel::detectCores())
options(mc.cores = 4)
rstan_options(auto_write = TRUE)

source("code/ctp_data_final.R")

source("code/Indiana_final.R")

source("code/Ohio_final.R")

source("code/BHM_final.R")

source("code/OtherStates_final.R")


