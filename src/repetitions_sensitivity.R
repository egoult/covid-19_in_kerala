##
## Check the sensitivity of the model to the number of repetitions
## Author: Elizabeth Goult
## Date: 11/01/2020
## Comments:


# clear work space
# clear graphics
rm(list=ls())
graphics.off()

# library
require(deSolve)
require(ggplot2)
require(tidyverse)
require(parallel)



cl<-makeCluster(detectCores()-1)

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

# functions

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

QModel<-function(time, X, pars){
    # X compartments 
    # S,E,I,R
    # SQ, EQ, IQ, RQ
    # Fitted kerala model
    

    with(as.list(c(X, pars)),{
        H<-sum(c(S,E,I,R))

        # w<- 0.000
        # omega <-1
        # omegaw<-1
        # specificity<- 1
        # sensitivity<-0.85
        

        if(time<(as.Date("2020-03-23") - as.Date(start_date))){
            total_delta<-round(rnorm(1, 46000, 2000))
            ninf<-rbinom(1, size=total_delta, prob=global[which((global_date-as.Date(start_date))==floor(time))]/7.8e+9) 

            deltas<-total_delta-ninf
            deltae<-runif(1,0, ninf)
            deltai<-runif(1,0, ninf  - deltae)
            deltar<-ninf  - deltae - deltai
       }else if(time<(as.Date("2020-05-16") - as.Date(start_date))){
            total_delta<-round(runif(1,0,920*2))
            ninf<-rbinom(1, size=total_delta, prob=global[which((global_date-as.Date(start_date))==floor(time))]/7.8e+9)
            
            deltas<-total_delta-ninf
            deltae<-runif(1,0, ninf)
            deltai<-runif(1,0, ninf  - deltae)
            deltar<-ninf  - deltae - deltai
        }else{
            total_delta<-travelin[which((travelin_date-as.Date(start_date))==floor(time))]
            ninf<-rbinom(1, size=total_delta, prob=global[which((global_date-as.Date(start_date))==floor(time))]/7.8e+9)
            
            deltas<-total_delta-ninf
            deltae<-runif(1,0, ninf)
            deltai<-runif(1,0, ninf  - deltae)
            deltar<-ninf  - deltae - deltai
            }

        if(time<(as.Date("2020-03-23") - as.Date(start_date))){
            lambda<-lambda1
        }else if (time<(as.Date("2020-04-20") - as.Date(start_date))){ 
            lambda<-lambda2
        }else if (time<(as.Date("2020-04-24") - as.Date(start_date))){ 
            lambda<-lambda3
        }else{
            lambda<-lambda3}

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

solveOde<-function(x, parameters){
    # parameters - named vector of parameters
    return(ode(y=State,Mod_times, func=QModel, parms=parameters, method="rk4"))
}

ConvertToArray<-function(x){
    # x - list of dataframes of the same size
    nlist<-length(x)

    as_array<-array(data=NA, dim=c(nrow(x[[1]]), ncol(x[[1]]), nlist), dimnames=list(NULL, colnames(x[[1]]), NULL))
    for(i in 1:nlist){as_array[,,i]<-x[[i]]}
    return(as_array)
}

ModelWrapper<-function(nrep){
    #par - vector of parameters
    # run the model
    x<-ConvertToArray(parLapply(cl=cl, X=1:nrep, fun=solveOde, parameters=parameters))

    # calculate mean values over runs
    meanx<-apply(x, c(1,2), mean)
    rm(x)

    #pull out max deaths, cumulative cases and cumulative hosp
    deaths<-floor(meanx[nrow(meanx),"Death"])
    cum_cases<-meanx[nrow(meanx), "New"]
    cum_hosp<-meanx[nrow(meanx), "New_hosp"]

    return(c(deaths, cum_cases, cum_hosp))}

DFWrapper<-function(par_df){
    # par_df - dataframe of parameters
    out<-vector(mode = "list", length = nrow(par_df))
    for(i in 1: nrow(par_df)){
        out[[i]]<-ModelWrapper(par=par_df[i,])
    }
    df<-as.data.frame(do.call(rbind, out))

    return(df)
}

# End of functions

# All model inputs
State<-c(S=33300000,E=0,I=0,R=0,SQ=0,EQ=0,IQ=0,RQ=0,Death = 0, New=0, New_hosp=0)
Mod_times<- seq(0.1,as.double(as.Date(end_date) - as.Date(start_date)), stepsize) 
parameters<-c(lambda1 = 1.12999975, lambda2=0.10999994, lambda3 = 1.09000037, 
            w= 0.000, omega = 1, p=0.19999992, r=0.07142859, 
            sigma=0.44220802, omegaw=1, d=0.0004800 , specificity = 1, sensitivity=0.85)

clusterExport(cl = cl, varlist=c("State", "Mod_times", "parameters", "ode", 
        "QModel", "start_date", "global", "global_date", "travelin", "travelin_date"))

nrep_values<-c(1, 3, 10, 30, 100, 300, 1000)
ntest<-30
#store results
res<-matrix(data=NA,ncol=5, nrow=ntest*length(nrep_values))
colnames(res)<-c("nrep", "test_no", "deaths", "cum_cases", "cum_hosp")
rowcounter<-0
tick<-Sys.time()
for (i in 1:length(nrep_values)){
    for(j in 1:ntest){
        #print(nrep_values[i])
        rowcounter<-rowcounter+1
        res[rowcounter,]<-c(nrep_values[i], j, ModelWrapper(nrep=nrep_values[i]))
    }
}

stopCluster(cl)
rm(cl)
tock<-Sys.time()

print(tock - tick)

res <-as.data.frame(res)

write.csv(res, file = "results/repetitions_sensitivity.csv")

# plot
plt_deaths<-ggplot(res, aes(x=nrep, y=deaths))+
    geom_point()+
    theme_classic()+
    labs(x= "Number of model runs", y="Cumulative deaths on 30th May")

plt_cum_cases<-ggplot(res, aes(x=nrep, y=cum_cases))+
    geom_point()+
    theme_classic()+
    labs(x= "Number of repititions", y=expression(Cumulative~cases~on~30^{th}~May~2020))

plt_cum_hosp<-ggplot(res, aes(x=nrep, y=cum_hosp))+
    geom_point()+
    theme_classic()+
    labs(x= "Number of repititions", y=expression(Cumulative~hospitalisations~on~30^{th}~May~2020))

# calculate "true" value
best_val<-res %>%
    filter(nrep == max(nrep_values)) %>%
    summarise(avr_death = mean(deaths), avr_cases = mean(cum_cases), avr_hosp = mean(cum_hosp))


plt_rel_death<-ggplot(res, aes(x=as.factor(nrep), y = (deaths/best_val$avr_death - 1)))+
    geom_point()+
    theme_classic()+
    labs(x= "Number of repititions", y="Relative difference between estimate and max repetitions estimate")

print(plt_rel_death)
dev.new()

plt_rel_hosp<-ggplot(res, aes(x=as.factor(nrep), y = (cum_hosp/best_val$avr_hosp - 1)))+
    geom_point()+
    theme_classic()+
    labs(x= "Number of repititions", y="Relative difference between estimate and max repetitions estimate")

print(plt_rel_hosp)


