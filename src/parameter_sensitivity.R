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
require(ggplot2)
require(tidyverse)
require(reshape2)
# require(parallel)
require(FME)

# cl<-makeCluster(detectCores()-1)

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
nrep<-1
#sobolsize<-20

# functions

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

QModel<-function(time, X, pars){
    # X compartments 
    # S,E,I,R
    # SQ, EQ, IQ, RQ
    # Fitted kerala model

    with(as.list(c(X, pars)),{
        H<-sum(c(S,E,I,R))

         w<- 0.000
        # omega <-1
        # omegaw<-1
        # specificity<- 1
        # sensitivity<-0.85
        

        if(time<(as.Date("2020-03-23") - as.Date(start_date))){
            total_delta<-46000#<-round(rnorm(1, 46000, 2000))
            #ninf<-rbinom(1, size=total_delta, prob=global[which((global_date-as.Date(start_date))==floor(time))]/7.8e+9) 

            # deltas<-total_delta-ninf
            # deltae<-ninf/3#runif(1,0, ninf)
            # deltai<-(ninf-deltae)/2#runif(1,0, ninf  - deltae)
            # deltar<-ninf  - deltae - deltai
       }else if(time<(as.Date("2020-05-16") - as.Date(start_date))){
            total_delta<-920#round(runif(1,0,920*2))
            #ninf<-rbinom(1, size=total_delta, prob=global[which((global_date-as.Date(start_date))==floor(time))]/7.8e+9)
            
            # deltas<-total_delta-ninf
            # deltae<-runif(1,0, ninf)
            # deltai<-runif(1,0, ninf  - deltae)
            # deltar<-ninf  - deltae - deltai
        }else{
            total_delta<-travelin[which((travelin_date-as.Date(start_date))==floor(time))]
            #ninf<-rbinom(1, size=total_delta, prob=global[which((global_date-as.Date(start_date))==floor(time))]/7.8e+9)
            
            # deltas<-total_delta-ninf
            # deltae<-runif(1,0, ninf)
            # deltai<-runif(1,0, ninf  - deltae)
            # deltar<-ninf  - deltae - deltai
            }
            glob<-global[which((global_date-as.Date(start_date))==floor(time))]/7.8e+9
            ninf<- total_delta*glob

            deltas<-total_delta-ninf
            deltae<-ninf/3#runif(1,0, ninf)
            deltai<-(ninf-deltae)/2#runif(1,0, ninf  - deltae)
            deltar<-ninf  - deltae - deltai
            
        if(time<(as.Date("2020-03-23") - as.Date(start_date))){
            lambda<-lambda1
        }else if (time<(as.Date("2020-04-20") - as.Date(start_date))){ 
            lambda<-lambda2
        }else if (time<(as.Date("2020-04-24") - as.Date(start_date))){ 
            lambda<-lambda3
        }else{
            lambda<-lambda3}

        #Non hospital compartments
        ds <- -lambda*S*I/H + w*R + omegaw*SQ + specificity*deltas-total_delta
        de <-  lambda*S*I/H - p*E + (1-sensitivity)*deltae 
        di <-  p*E - r*I - sigma*I + (1-sensitivity)*deltai
        dr <-  r*I + omegaw*RQ - w*R +specificity* deltar

        # hospital
        # work out number of deaths
        
        dailydeaths<-IQ*d#rbinom(1,size=round(IQ), prob=d)
        #dailydeaths<-rbinom(1,size=round(IQ), prob=d)
        
        dsq <- (1-specificity)*deltas- omegaw*SQ 
        deq <- sensitivity*deltae- p*EQ 
        diq <- p*EQ + sigma*I + sensitivity*deltai -r*IQ - dailydeaths
        drq <- r*IQ - omegaw*RQ +(1-specificity)*deltar

        ddeath<-dailydeaths

        dnew_cases<-sum(c(p*E + (1-sensitivity)*deltai, p*EQ + sensitivity*deltai))
        dnew_hosp<-sum(c((1-specificity)*deltas,sensitivity*deltae, sigma*I + sensitivity*deltai,(1-specificity)*deltar))
    return(list(c(ds,de,di,dr,dsq,deq,diq,drq, ddeath, dnew_cases, dnew_hosp)))})}

# solveOde<-function(x, parameters){
#     # parameters - named vector of parameters
#     return(ode(y=State,Mod_times, func=QModel, parms=parameters, method="rk4"))
# }
solveOde<-function(parameters){
    # parameters - named vector of parameters
    return(ode(y=State,Mod_times, func=QModel, parms=parameters, method="rk4"))
}



ConvertToArray<-function(x){
    # x - list of dataframes of the same size
    nlist<-length(x)

    as_array<-array(data=NA, dim=c(nrow(x[[1]]), ncol(x[[1]]), nlist), dimnames=list(NULL, colnames(x[[1]]), NULL))
    for(i in 1:nlist){as_array[,,i]<-x[[i]]}
    rm(x)
    return(as_array)
}

ModelWrapper<-function(parms){
    #par - vector of parameters
    # run the model
    
    meanx<-solveOde(parms)
    # x<-ConvertToArray(parLapply(cl=cl, X=1:nrep, fun=solveOde, parameters=parms))

    # # calculate mean values over runs
    # meanx<-apply(x, c(1,2), mean)
    # rm(x)
    
    #pull out time, deaths, cumulative cases and cumulative hosp
    time<-meanx[is.wholenumber(meanx[,"time"]),"time"]
    deaths<-meanx[is.wholenumber(meanx[,"time"]),"Death"]
    cum_cases<-meanx[is.wholenumber(meanx[,"time"]), "New"]
    cum_hosp<-meanx[is.wholenumber(meanx[,"time"]), "New_hosp"]

    return(as.data.frame(cbind(time, deaths, cum_cases, cum_hosp)))}

# DFWrapper<-function(par_df){
#     # par_df - dataframe of parameters
#     out<-vector(mode = "list", length = nrow(par_df))
#     for(i in 1: nrow(par_df)){
#         out[[i]]<-ModelWrapper(par=par_df[i,])
#     }
#     df<-as.data.frame(do.call(rbind, out))

#     return(df)
# }

# End of functions

# All model inputs
State<-c(S=33300000,E=0,I=0,R=0,SQ=0,EQ=0,IQ=0,RQ=0,Death = 0, New=0, New_hosp=0)
Mod_times<- seq(0.1,as.double(as.Date(end_date) - as.Date(start_date)), stepsize) 
parameters<-list(lambda1 = 1.12999975, lambda2=0.10999994, lambda3 = 1.09000037, 
              p=0.19999992, r=0.07142859, 
            sigma=0.44220802, omegaw=1, d=0.0004800 , specificity = 1, sensitivity=0.85)

parameters_min <- c(lambda1 = 0, lambda2= 0, lambda3 = 0, 
            w= 0.000, omega = 1, p=0, r=0, 
            sigma=0, omegaw=1, d=0, specificity = 0, sensitivity=0)
parameters_max<-c(lambda1 = 2, lambda2=2, lambda3 = 2, 
            w= 0.000, omega = 1, p=1, r=2, 
            sigma=1, omegaw=1, d=1 , specificity = 1, sensitivity=1)
#stop()
# clusterExport(cl=cl, varlist=c("State", "Mod_times", "parameters", "ode", 
#         "QModel", "start_date", "global", "global_date", "travelin", "travelin_date", "ModelWrapper"))


factors <- c("lambda1", "lambda2", "lambda3", 
            "w", "omega", "p", "r", 
            "sigma", "omegaw", "d" , 
            "specificity", "sensitivity")

sa<-sensFun(func=ModelWrapper, parms=parameters, sensvar = c("deaths", "cum_cases"),#, "cum_hosp"), 
            varscale = NULL, parscale = NULL, tiny = 1e-8, map = 1)

write.csv(sa, "results/parameter_sa_det.csv")
dev.new()
pairs(sa)
dev.new()
plot(sa)

samelt<-melt(sa,id.vars=c("x", "var"))

dev.new()
plt<-ggplot(samelt, aes(x=x,y=value)) +
    geom_line() +
    facet_grid(var~variable, scales="free_y")
print(plt)

sa_deaths<-sa %>%
        filter(var=="deaths")
write.csv(summary(sa_deaths), "results/sa_deaths.csv")
write.csv(sa_deaths[nrow(sa_deaths),], "results/sa_deaths_max.csv")

sa_cum_cases<-sa %>%
    filter(var=="cum_cases")
write.csv(summary(sa_cum_cases), "results/sa_cum_cases.csv")
write.csv(sa_cum_cases[nrow(sa_cum_cases),], "results/sa_cum_cases_max.csv")

samelt_deaths<-samelt %>%
    filter(var=="deaths")

plt_deaths<-ggplot(samelt_deaths, aes(x=x, y=value))+
     geom_line(aes(color = variable, linetype = variable)) +
     theme_classic()+
     labs(x="Days since first observed case", y="Sensitivity to parameter")


samelt_cum_cases<-samelt %>%
    filter(var=="cum_cases")

plt_cum_cases<-ggplot(samelt_cum_cases, aes(x=x, y=value))+
     geom_line(aes(color = variable, linetype = variable)) +
     theme_classic()+
     labs(x="Days since first observed case", y="Sensitivity to parameter")
