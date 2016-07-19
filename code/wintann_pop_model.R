####################################################################
##  wintann_pop_model.R: simulates long term population dynamics  ##
##  of species in a Sonoran Desert winter annual plant community. ##
##  The model is a generic annual plant model with seedbank.      ##
####################################################################

# Clear the workspace
rm(list=ls(all.names = TRUE))

# Load some libraries
library(plyr)
library(reshape2)
library(ggthemes)
library(ggplot2)

do_species <- "erla"

####
####  ANNUAL PLANT POPULATION MODEL
####
project_popmod <- function(seed_surv, seeds_per_seedling, germination_frac)
{
  pgr <- seed_surv*(1-germination_frac) + seeds_per_seedling*germination_frac
}



####
####  READ IN ANNUAL PLANT DATA
####
species_year_data <- read.csv("../data/species_x_year.csv")
single_spp_data <- subset(species_year_data, species==do_species)



####
####  MODEL LONG-TERM PER CAPITA GROWTH RATE
####
single_spp_data <- subset(single_spp_data, germ.fraction!="." & lxbx!=".")
seeds_per_seedling_vec <- as.numeric(as.character(single_spp_data$lxbx))
seed_surv <- 0
germination_frac_vec <- as.numeric(as.character(single_spp_data$germ.fraction))
iters <- 100
pgr <- numeric(iters)
for(t in 1:iters){
  randyr <- sample(c(1:length(seeds_per_seedling_vec)), 1)
  pgr[t] <- project_popmod(seed_surv = seed_surv, 
                        seeds_per_seedling = seeds_per_seedling_vec[randyr],
                        germination_frac = germination_frac_vec[randyr])
}

plot(pgr, type="l")
