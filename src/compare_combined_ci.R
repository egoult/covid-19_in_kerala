##
## Plot to compare the combined resamplingvs ci resampling results.
## Author: elg
## Comments:

# clear workspace
# clear graphics
rm(list = ls())
graphics.off()

# packages
require(tidyverse)
require(ggplot2)

# data
start_date<- "2020-01-30" #"2020-03-08"
end_date  <- "2020-05-30" #"2020-05-30"

kerala<-read.csv("RAW_DATA/kerala_covid19_20200712.csv")
obs_date<-as.Date(kerala$Date, "%d/%m/%y")
cases<-kerala$Current_cases[which(obs_date==start_date):which(obs_date==end_date)] 
deaths<-kerala$Cumulative_deaths[which(obs_date==start_date):which(obs_date==end_date)]
obs<-data.frame(Date = as.Date(start_date) + seq(0, length(deaths)-1), cases = cases, deaths = deaths) %>%
    mutate(type = "Observations")

# simulations
ci_sim<-read.csv("results/CI/hospitalised_cases_uncertainty.csv") %>%
    mutate(type = "ci")
comb_sim<-read.csv("results/combined/hospitalised_cases_uncertainty.csv") %>%
    mutate(type = "combined")

sim<-bind_rows(ci_sim, comb_sim) %>%
    mutate(Date = as.Date(start_date)+Time) %>%
    select(!Time)




plt<-ggplot(obs, aes(x = Date, y = cases, group = type))+
    geom_point()+
    geom_line(data = sim, aes(x=Date, y = Mean, group = type, color = type))+
    geom_line(data = sim, aes(x=Date, y = Min, group = type, color = type))+
    geom_line(data = sim, aes(x=Date, y = Max, group = type, color = type))



