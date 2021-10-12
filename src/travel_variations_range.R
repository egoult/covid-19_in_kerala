##
## Produces the hospitalisation model with reduced testing
## And returns a csv with output
## Author: Elizabeth Goult
## Date: 08/05/2021
## Comments: Models
##          Original quarantine model
##          Reduced testing rate


# clear work space
# clear graphics
rm(list=ls())
graphics.off()

# library
require(deSolve)
require(parallel)
require(miscset)
require(tidyverse)
require(ggplot2)
require(viridis)
require(patchwork)

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

reporting_delay<-7
stepsize<-0.1
nrep<-100

#Functions
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

QModel<-function(time, X, pars){
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


ReducedQModel<-function(time, X, pars){
    # X compartments 
    # S,E,I,R
    # SQ, EQ, IQ, RQ
    # Testing rate reduced from 100% to a 

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
            lambda<-lambda4}

        #Non hospitalised compartments
        ds <- -lambda*S*I/H + w*R + omega*SQ + specificity*a*deltas+(1-a)*deltas-total_delta
        de <-  lambda*S*I/H - p*E + (1-sensitivity)*a*deltae + (1-a)*deltae 
        di <-  p*E - r*I - sigma*I + (1-sensitivity)*a*deltai + (1-a)*deltai
        dr <-  r*I + omegaw*RQ - w*R + a*specificity*deltar+(1-a)*deltar

        # hospitalised
        # work out number of deaths
        
        dailydeaths<-rbinom(1,size=round(IQ), prob=d)
        
        dsq <- (1-specificity)*a*deltas- omega*SQ 
        deq <- sensitivity*a*deltae- p*EQ 
        diq <- p*EQ + sigma*I + sensitivity*a*deltai -r*IQ - dailydeaths
        drq <- r*IQ - omegaw*RQ + a*(1-specificity)*deltar

        ddeath<-dailydeaths
        dnew_cases<-sum(c(p*E + (1-sensitivity)*a*deltai + (1-a)*deltai, p*EQ + sensitivity*a*deltai))
        dnew_hosp<-sum(c((1-specificity)*a*deltas,sensitivity*a*deltae, sigma*I + sensitivity*a*deltai,(1-specificity)*a*deltar))
    return(list(c(ds,de,di,dr,dsq,deq,diq,drq, ddeath, dnew_cases, dnew_hosp)))})}

    

    # NRMSD<-function(mod, obs){
    #     sqdiff<-(mod-obs)^2
    #     rootmn<-sqrt(mean(sqdiff))
    #     norm<-rootmn/mean(obs)
    #     return(norm)}


    NewCases<-function(output){
        I<-output[,4]+output[,8]
        R<-output[,5]+output[,9]
        deaths<-output[,10]

        It<-c(I,NA)
        Itm1<-c(NA,I)

        Rt<-c(R,NA)
        Rtm1<-c(NA,R)

        deathst<-c(deaths, NA)
        deathstm1<-c(NA, deaths)

        newI<-It-Itm1 - (Rt - Rtm1) - (deathst - deathstm1)
        return(newI)}

Reporting<-function(solve, solveCI, hosp, total){
    Active_total_cases<-rowSums(solve[,total])
    Active_total_cases_CI<-rowSums(solveCI[,total])

    Cumulative_total_cases<-solve[,"New"]
    Cumulative_total_cases_CI<-solveCI[,"New"]
    print(paste("Total cases", max(Active_total_cases), max(Cumulative_total_cases)))
    print(paste("Total cases CI +", max(Active_total_cases+Active_total_cases_CI), max(Cumulative_total_cases+Cumulative_total_cases_CI)))
    print(paste("Total cases CI -", max(Active_total_cases-Active_total_cases_CI), max(Cumulative_total_cases-Cumulative_total_cases_CI)))

    Active_hospital_cases<-rowSums(solve[,hosp])
    Active_hospital_cases_CI<-rowSums(solveCI[,hosp])
    
    Cumulative_hospital_cases<-solve[,"New_hosp"]
    Cumulative_hospital_cases_CI<-solveCI[,"New_hosp"]

    shifted_time<-solve[,"time"]+reporting_delay
    shifted<-which(shifted_time<=max(solve[,"time"]))
    print(paste("Reporting shifted hospitalised cases", max(Active_hospital_cases[shifted]), max(Cumulative_hospital_cases[shifted])))
    print(paste("Reporting shifted hospitalised cases CI +", max(Active_hospital_cases[shifted]+Active_hospital_cases_CI[shifted]), max(Cumulative_hospital_cases[shifted]+Cumulative_hospital_cases_CI[shifted])))
    print(paste("Reporting shifted hospitalised cases CI -", max(Active_hospital_cases[shifted]-Active_hospital_cases_CI[shifted]), max(Cumulative_hospital_cases[shifted]-Cumulative_hospital_cases_CI[shifted])))

    Deaths<-solve[,"Death"]
    Deaths_CI<-solveCI[,"Death"]
    print(paste0("Deaths ",max(Deaths)))
    print(paste0("Deaths + CI ",max(Deaths+Deaths_CI)))
    print(paste0("Deaths - CI ",max(Deaths-Deaths_CI)))    

    print(paste("mortality rates", max(Deaths)/max(Cumulative_hospital_cases[shifted]), max(Deaths)/max(Cumulative_total_cases) ))
    print(paste("mortality rates", max(Deaths+Deaths_CI)/max(Cumulative_hospital_cases[shifted]-Cumulative_hospital_cases_CI[shifted]), max(Deaths+Deaths_CI)/max(Cumulative_total_cases-Cumulative_total_cases_CI) ))
    print(paste("mortality rates", max(Deaths-Deaths_CI)/max(Cumulative_hospital_cases[shifted]+Cumulative_hospital_cases_CI[shifted]), max(Deaths-Deaths_CI)/max(Cumulative_total_cases+Cumulative_total_cases_CI) ))
    df<-data.frame(object=c("max_cases", "cum_cases", "max_hosp", "cum_hosp", "deaths", "cfr", "ifr"),
                   value=c(max(Active_total_cases), max(Cumulative_total_cases),
                           max(Active_hospital_cases[shifted]), max(Cumulative_hospital_cases[shifted]),
                           max(Deaths), max(Deaths)/max(Cumulative_hospital_cases[shifted]), max(Deaths)/max(Cumulative_total_cases)),
                   low_95=c(max(Active_total_cases-Active_total_cases_CI), max(Cumulative_total_cases-Cumulative_total_cases_CI),
                            max(Active_hospital_cases[shifted]-Active_hospital_cases_CI[shifted]), max(Cumulative_hospital_cases[shifted]-Cumulative_hospital_cases_CI[shifted]),
                            max(Deaths-Deaths_CI), max(Deaths-Deaths_CI)/max(Cumulative_hospital_cases[shifted]+Cumulative_hospital_cases_CI[shifted]), max(Deaths-Deaths_CI)/max(Cumulative_total_cases+Cumulative_total_cases_CI)),
                    up_95=c(max(Active_total_cases+Active_total_cases_CI), max(Cumulative_total_cases+Cumulative_total_cases_CI),
                            max(Active_hospital_cases[shifted]+Active_hospital_cases_CI[shifted]), max(Cumulative_hospital_cases[shifted]+Cumulative_hospital_cases_CI[shifted]),
                            max(Deaths+Deaths_CI), max(Deaths+Deaths_CI)/max(Cumulative_hospital_cases[shifted]-Cumulative_hospital_cases_CI[shifted]), max(Deaths+Deaths_CI)/max(Cumulative_total_cases-Cumulative_total_cases_CI)))
    rm(solve)
    return(df)
}

ConvertToArray<-function(x){
    # x - list of dataframes of the same size
    nlist<-length(x)

    as_array<-array(data=NA, dim=c(nrow(x[[1]]), ncol(x[[1]]), nlist), dimnames=list(NULL, colnames(x[[1]]), NULL))
    for(i in 1:nlist){as_array[,,i]<-x[[i]]}
    return(as_array)
}
geq1<-function(x){
    if(x>1){return(x)}
    return(1)}
GEQvector<-function(x){
    sapply(x, geq1)
}

TestingRateSolve<-function(testing_parameters){
    ## Solve the model for a given testing rate contained in parameters
    ## Inputs:
    ## rate: numeric, testing rate
    ## parameters: named list of parameters


    clusterExport(cl = cl, varlist=c("testing_parameters","ReducedQModel"))

    testing_solve_array<-ConvertToArray(parLapply(cl=cl, X=1:nrep, fun = function(x) ode(y=State, Mod_times, func=ReducedQModel, parms=testing_parameters, method="rk4")))
    # testing_solve<-apply(testing_solve_array, c(1,2), mean)
    # testing_solve_CI<-apply(testing_solve_array, c(1,2), FUN=function(x)confint(x, parm=qnorm)) # half of CI
    # testing_solve_death<-floor(testing_solve[,"Death"])[which(is.wholenumber(Mod_times))]

    return(testing_solve_array)
}


# end functions

# All model inputs
State<-c(S=33300000,E=0,I=0,R=0,SQ=0,EQ=0,IQ=0,RQ=0,Death = 0, New=0, New_hosp=0)
Mod_times<- seq(0,as.double(as.Date(end_date) - as.Date(start_date)), stepsize) 
parameters<-c(lambda1 = 1.2268909911, lambda2 = 0.1657483443, lambda3=1.2470045370, lambda4=1.2470045370,  sigma= 0.4971470720, d= 0.0004840214, p= 0.2293377863,r= 0.1042903831, cases_report= 7.0000000000 )

clusterExport(cl = cl, varlist=c("State", "Mod_times", "parameters", "ode", 
                                "QModel","start_date", "global", "global_date", "travelin", 
                                "travelin_date"))



# quarantine_solve_array<-ConvertToArray(parLapply(cl=cl, X=1:nrep, fun = function(x) ode(y=State, Mod_times, func=QModel, parms=parameters, method="rk4")))

# quarantine_solve<-apply(quarantine_solve_array, c(1,2), mean)
# quarantine_solve_CI<-apply(quarantine_solve_array, c(1,2), FUN=function(x)confint(x, parm=qnorm)) # half of CI
# quarantine_solve_death<-floor(quarantine_solve[,"Death"])[which(is.wholenumber(Mod_times))]


#Reduced testing rate 10%
# clusterExport(cl = cl, varlist=c("testing_parameters","ReducedQModel"))
# #testing_solve<-replicate(n=nrep, ode(y=State,Mod_times,func=ReducedQModel, parms=testing_parameters, method="rk4") )
# testing_solve_array<-ConvertToArray(parLapply(cl=cl, X=1:nrep, fun = function(x) ode(y=State, Mod_times, func=ReducedQModel, parms=testing_parameters, method="rk4")))
# testing_solve<-apply(testing_solve_array, c(1,2), mean)
# testing_solve_CI<-apply(testing_solve_array, c(1,2), FUN=function(x)confint(x, parm=qnorm)) # half of CI
# testing_solve_death<-floor(testing_solve[,"Death"])[which(is.wholenumber(Mod_times))]

testing_parameters<-c(parameters,a=0)
t0<- TestingRateSolve(testing_parameters) %>%
    as.data.frame.table() %>%
    pivot_wider( names_from = "Var2", values_from = "Freq") %>%
    select(!Var1) %>%
    mutate(Rate = 0)


testing_parameters<-c(parameters,a=0.1)
t10<- TestingRateSolve(testing_parameters) %>%
    as.data.frame.table() %>%
    pivot_wider( names_from = "Var2", values_from = "Freq") %>%
    select(!Var1) %>%
    mutate(Rate = 10)

testing_parameters<-c(parameters,a=0.2)
t20<- TestingRateSolve(testing_parameters) %>%
    as.data.frame.table() %>%
    pivot_wider( names_from = "Var2", values_from = "Freq") %>%
    select(!Var1) %>%
    mutate(Rate = 20)

testing_parameters<-c(parameters,a=0.3)
t30<- TestingRateSolve(testing_parameters) %>%
    as.data.frame.table() %>%
    pivot_wider( names_from = "Var2", values_from = "Freq") %>%
    select(!Var1) %>%
    mutate(Rate = 30)


testing_parameters<-c(parameters,a=0.4)
t40<- TestingRateSolve(testing_parameters) %>%
    as.data.frame.table() %>%
    pivot_wider( names_from = "Var2", values_from = "Freq") %>%
    select(!Var1) %>%
    mutate(Rate = 40)


testing_parameters<-c(parameters,a=0.5)
t50<- TestingRateSolve(testing_parameters) %>%
    as.data.frame.table() %>%
    pivot_wider( names_from = "Var2", values_from = "Freq") %>%
    select(!Var1) %>%
    mutate(Rate = 50)


testing_parameters<-c(parameters,a=0.6)
t60<- TestingRateSolve(testing_parameters) %>%
    as.data.frame.table() %>%
    pivot_wider( names_from = "Var2", values_from = "Freq") %>%
    select(!Var1) %>%
    mutate(Rate = 60)

testing_parameters<-c(parameters,a=0.7)
t70<- TestingRateSolve(testing_parameters) %>%
    as.data.frame.table() %>%
    pivot_wider( names_from = "Var2", values_from = "Freq") %>%
    select(!Var1) %>%
    mutate(Rate = 70)

testing_parameters<-c(parameters,a=0.8)
t80<- TestingRateSolve(testing_parameters) %>%
    as.data.frame.table() %>%
    pivot_wider( names_from = "Var2", values_from = "Freq") %>%
    select(!Var1) %>%
    mutate(Rate = 80)


testing_parameters<-c(parameters,a=0.9)
t90<- TestingRateSolve(testing_parameters) %>%
    as.data.frame.table() %>%
    pivot_wider( names_from = "Var2", values_from = "Freq") %>%
    select(!Var1) %>%
    mutate(Rate = 90)

testing_parameters<-c(parameters,a=1)
t100<- TestingRateSolve(testing_parameters) %>%
    as.data.frame.table() %>%
    pivot_wider( names_from = "Var2", values_from = "Freq") %>%
    select(!Var1) %>%
    mutate(Rate = 100)



testing<-rbind(t0,t10,t20,t30,t40,t50,t60,t70,t80,t90, t100)

stopCluster(cl)
rm(cl)

testing_mean<-testing %>%
    group_by(Rate, time) %>%
    summarise(Deaths = round(mean(Death)), 
            reported_cases = mean(SQ+EQ+IQ+RQ), 
            actual_cases = mean(I+IQ), 
            Death_CI = confint(Death, parm=qnorm)[1],
            reported_cases_CI = confint(SQ+EQ+IQ+RQ, parm=qnorm)[1], 
            actual_cases_CI = confint(I+IQ, parm=qnorm)[1] ) %>%
    mutate(Rate = 0.01*Rate) 


plt.reported<-ggplot(testing_mean, aes(x = time, y = reported_cases, group = Rate, color=Rate))+
   geom_line()+
   geom_ribbon(aes(ymin = reported_cases-reported_cases_CI, ymax = reported_cases+reported_cases_CI, fill = Rate), alpha = 0.7)+
   theme_classic()+
   scale_fill_manual(values = viridis(11))+
   labs(x = "Days since initial infection", y = "Reported cases")


plt.actual<-ggplot(testing_mean, aes(x = time, y = actual_cases, group = Rate, color=Rate))+
   geom_line()+
   geom_ribbon(aes(ymin = actual_cases-actual_cases_CI, ymax = actual_cases+actual_cases_CI, color = Rate, fill = Rate), alpha = 0.7)+
   theme_classic()+
   labs(x ="Days since initial infection", y = "Total cases")#+
   #scale_fill_manual(values = wes_palette("GrandBudapest1", n = 10, type = "continuous"))

plt.death<-ggplot(testing_mean, aes(x = time, y = Deaths, group = Rate, color=Rate))+
   geom_line()+
   geom_ribbon(aes(ymin = round(Deaths-Death_CI), ymax = round(Deaths+Death_CI), color = Rate, fill = Rate), alpha = 0.7)+
   theme_classic()+
   labs(x ="Days since initial infection", y = "Deaths")


print(plt.actual + plt.death)