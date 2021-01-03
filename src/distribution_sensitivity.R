##
## Check the impact of the choice of distribution for travel in and mortality
## Author: Elizabeth Goult
## Date : 03/01/2020
## comments:

# clear work space
# clear graphics
rm(list=ls())
graphics.off()

# library
require(deSolve)
require(ggplot2)
require(tidyverse)
require(reshape2)

# read in global timeseries of number of cases of covid19
global<-read.csv("RAW_DATA/global_timeseries_20200712.csv")
global_date<-as.Date(global$Date, "%d/%m/%y")
global<-global$Global[!is.na(global$Global)]

# travel into the state
travelin<-read.csv("RAW_DATA/travel_into_kerala_20200712.csv")
travelin_date<-as.Date(travelin$Date, "%d/%m/%y")
travelin<-travelin$Daily_travel

# kerala data for comparison
kerala<-read.csv("RAW_DATA/kerala_covid19_20200712.csv")
kcases<-kerala$Current_cases
kdeaths<-kerala$Cumulative_deaths


start_date<- "2020-01-30" #"2020-03-08"
end_date  <- "2020-05-30"

reporting_delay<-4
stepsize<-0.1
nrep<-3

dists<-c("uniform", "normal")

# Functions
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

QModel<-function(time, X, pars, eir_fun){
    # X compartments 
    # S,E,I,R
    # SQ, EQ, IQ, RQ
    # Fitted kerala model
    

    with(as.list(c(X, pars)),{
        H<-sum(c(S,E,I,R))

        w<- 0.000
        omega <-1
        omegaw<-1
        specificity<- 1
        sensitivity<-0.85
        

        if(time<(as.Date("2020-03-23") - as.Date(start_date))){
            total_delta<-round(rnorm(1, 46000, 2000))
            ninf<-rbinom(1, size=total_delta, prob=global[which((global_date-as.Date(start_date))==floor(time))]/7.8e+9) 
       }else if(time<(as.Date("2020-05-16") - as.Date(start_date))){
            total_delta<-round(runif(1,0,920*2))
            ninf<-rbinom(1, size=total_delta, prob=global[which((global_date-as.Date(start_date))==floor(time))]/7.8e+9)
        }else{
            total_delta<-travelin[which((travelin_date-as.Date(start_date))==floor(time))]
            ninf<-rbinom(1, size=total_delta, prob=global[which((global_date-as.Date(start_date))==floor(time))]/7.8e+9)            
            }
        deltas<-total_delta-ninf

        if(!any(eir_fun %in% c("uniform", "normal"))){
            stop("eir_fun should be one of, 'uniform', 'normal'")}

        if(eir_fun == "uniform"){
            deltae<-runif(n=1,0, ninf)
            deltai<-runif(n=1,0, ninf - deltae)
            deltar<-ninf  - deltae - deltai}
        if(eir_fun == "normal"){
            deltae<-rnorm(n=1, mean=ninf/3, sd=ninf/10)
            deltai<-rnorm(n=1, mean=(ninf-deltae)/2, sd=ninf/10)
            deltar<-ninf  - deltae - deltai}
        
        if(time<(as.Date("2020-03-23") - as.Date(start_date))){
            lambda<-lambda1
        }else if (time<(as.Date("2020-04-20") - as.Date(start_date))){ 
            lambda<-lambda2
        }else if (time<(as.Date("2020-04-24") - as.Date(start_date))){ 
            lambda<-lambda3
        }else{
            lambda<-lambda4}

        #Non hospital compartments
        ds <- -lambda*S*I/H + w*R + omega*SQ + specificity*deltas-total_delta
        de <-  lambda*S*I/H - p*E + (1-sensitivity)*deltae 
        di <-  p*E - r*I - sigma*I + (1-sensitivity)*deltai
        dr <-  r*I + omegaw*RQ - w*R +specificity* deltar

        # hospital
        # work out number of deaths
        
        dailydeaths<-rbinom(1,size=round(IQ), prob=d)
        
        dsq <- (1-specificity)*deltas- omega*SQ 
        deq <- sensitivity*deltae- p*EQ 
        diq <- p*EQ + sigma*I + sensitivity*deltai -r*IQ - dailydeaths
        drq <- r*IQ - omegaw*RQ +(1-specificity)*deltar

        ddeath<-dailydeaths

        dnew_cases<-sum(c(p*E + (1-sensitivity)*deltai, p*EQ + sensitivity*deltai))
        dnew_hosp<-sum(c((1-specificity)*deltas,sensitivity*deltae, sigma*I + sensitivity*deltai,(1-specificity)*deltar))
    return(list(c(ds,de,di,dr,dsq,deq,diq,drq, ddeath, dnew_cases, dnew_hosp)))})}

# All model inputs
State<-c(S=33300000,E=0,I=0,R=0,SQ=0,EQ=0,IQ=0,RQ=0,Death = 0, New=0, New_hosp=0)
Mod_times<- seq(0,as.double(as.Date(end_date) - as.Date(start_date)), stepsize) 
parameters<-c(lambda1 = 1.12999975, lambda2=0.10999994, lambda3 = 1.09000037, lambda4 = 1.09000037, 
            w= 0.000, omega = 1, p=0.19999992, r=0.07142859, 
            sigma=0.44220802, omegaw=1, d=0.0004800 , specificity = 1, sensitivity=0.85)

# run for all distributions
dist_runs<-list()
for (i in 1:length(dists)){
    dist_runs[[i]]<-replicate(n=nrep, ode(y=State,Mod_times, func=QModel, parms=parameters, method="rk4", eir_fun=dists[i]) )
}

# plot results
# melt into one df



test<-melt(dist_runs[[i]], value.name = 'Score', id.vars="time") 