##
## Check the sensitivity of the model to parameter values
## Author: Elizabeth Goult
## Date: 05/01/2020
## Comments:


# clear work space
# clear graphics
rm(list=ls())
graphics.off()

# library
require(deSolve)
#require(ggplot2)
# require(tidyverse)
# require(reshape2)
require(parallel)
# require(ODEsensitivity)
require(pse)

cl<-makeCluster(1)#(detectCores()-1)

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
nrep<-20
sobolsize<-20

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

ModelWrapper<-function(par){
    #par - vector of parameters
    # run the model
    x<-ConvertToArray(parLapply(cl=cl, X=1:nrep, fun=solveOde, parameters=par))

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
tick<-Sys.time()
a<-ModelWrapper(par=parameters)
tock<-Sys.time()

print(tock - tick)



parameters_min <- c(lambda1 = 0, lambda2= 0, lambda3 = 0, 
            w= 0.000, omega = 1, p=0, r=0, 
            sigma=0, omegaw=1, d=0, specificity = 0, sensitivity=0)
parameters_max<-c(lambda1 = 2, lambda2=2, lambda3 = 2, 
            w= 0.000, omega = 1, p=1, r=2, 
            sigma=1, omegaw=1, d=1 , specificity = 1, sensitivity=1)

# clusterExport( varlist=c("State", "Mod_times", "parameters", "ode", 
#         "QModel", "start_date", "global", "global_date", "travelin", "travelin_date"))

# ODEsobol(mod =QModel, pars=parameters, state_init=State, times = Mod_times, n = sobolsize,
#     rfuncs = "runif", rargs = paste0("min = ", parameters_min,", max = ", parameters_min),
#     sobol_method = "Martinez", ode_method = "lsoda", 
#     parallel_eval = TRUE, parallel_eval_ncores = detectCores() - 1)


factors <- c("lambda1", "lambda2", "lambda3", 
            "w", "omega", "p", "r", 
            "sigma", "omegaw", "d" , 
            "specificity", "sensitivity")
quantiles <- c("qunif", "qunif", "qunif", 
            "qunif", "qunif", "qunif", "qunif", 
            "qunif", "qunif", "qunif" , 
            "qunif", "qunif")
q.args <- list(list(min=parameters_min[1], max=parameters_max[1]), 
            list(min=parameters_min[2], max=parameters_max[2]), 
            list(min=parameters_min[3], max=parameters_max[3]), 
            list(min=parameters_min[4], max=parameters_max[4]), 
            list(min=parameters_min[5], max=parameters_max[5]), 
            list(min=parameters_min[6], max=parameters_max[6]), 
            list(min=parameters_min[7], max=parameters_max[7]), 
            list(min=parameters_min[8], max=parameters_max[8]), 
            list(min=parameters_min[9], max=parameters_max[9]), 
            list(min=parameters_min[10], max=parameters_max[10]), 
            list(min=parameters_min[11], max=parameters_max[11]),
            list(min=parameters_min[12], max=parameters_max[12]))

myLHS <- LHS(DFWrapper, factors, N=20, quantiles, q.args, nboot=0, method="random")
print(plotecdf(myLHS, stack=TRUE))

print(plotscatter(myLHS))

# stopCluster(cl)
# rm(cl)