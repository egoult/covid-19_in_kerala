##
## Given a .rds of  an mcmc chain, runs the uncertainty
##
##
##


# clear workspace
rm(list = ls())
graphics.off()

# libraries
require(FME)
require(tidyverse)
# require(ggplot)

# fit set up
# read in global timeseries of number of cases of covid19
global<-read.csv("RAW_DATA/global_timeseries_20200712.csv")
global_date<-as.Date(global$Date, "%d/%m/%y")
global<-global$Global[!is.na(global$Global)]

# travel into the state
travelin<-read.csv("RAW_DATA/travel_into_kerala_20200712.csv")
travelin_date<-as.Date(travelin$Date, "%d/%m/%y")
travelin<-travelin$Daily_travel

start_date<- "2020-01-30" #"2020-03-08"
end_date  <- "2020-05-30" #"2020-05-30"
stepsize<-0.1
ndigits<-1
nrep<-1
nmcmc<-75000

nbin<-0

#Functions
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

QModel<-function(time, X, pars){
    # X compartments 
    # S,E,I,R
    # SQ, EQ, IQ, RQ
    # parameters

    with(as.list(c(X, pars)),{
        H<-sum(c(S,E,I,R))

        w<- 0.000
        omega <-1
        omegaw<-1
        specificity<- 1
        sensitivity<-0.85
        p<-0.2
        r<-1/14

        
        if(time<(as.Date("2020-03-24") - as.Date(start_date))){
            total_delta<-46000#round(rnorm(1, 46000, 2000))

        }else if(time<(as.Date("2020-05-16") - as.Date(start_date))){
            total_delta<-920#round(runif(1,0,920*2))
            
        }else{
            total_delta<-travelin[which((travelin_date-as.Date(start_date))==floor(time))]
            
            }

        global_rate <- global[which((global_date-as.Date(start_date))==floor(time))]/7.8e+9
        ninf<-total_delta * global_rate
        deltas<-total_delta-ninf
        deltae<-ninf/3#runif(1,0, ninf)
        deltai<-ninf/3#runif(1,0, ninf  - deltae)
        deltar<-ninf  - deltae - deltai

        if(time<(as.Date("2020-03-24") - as.Date(start_date))){
            lambda<-lambda1
        }else if (time<(as.Date("2020-04-20") - as.Date(start_date))){ 
            lambda<-lambda2
        }else if (time<(as.Date("2020-04-24") - as.Date(start_date))){ 
            lambda<-lambda3
        }else{
            lambda<-lambda3}
        
        
        #Non hospitalised compartments
        ds <- -lambda*S*I/H + w*R + omega*SQ + specificity*deltas-total_delta
        de <-  lambda*S*I/H - p*E + (1-sensitivity)*deltae 
        di <-  p*E - r*I - sigma*I + (1-sensitivity)*deltai
        dr <-  r*I + omegaw*RQ - w*R + specificity*deltar

        # hospitalised
        # work out number of deaths
        
        dailydeaths<-d*IQ #rbinom(1,size=round(IQ), prob=d)
        
        dsq <- (1-specificity)*deltas- omega*SQ 
        deq <- sensitivity*deltae- p*EQ 
        diq <- p*EQ + sigma*I + sensitivity*deltai -r*IQ - dailydeaths
        drq <- r*IQ - omegaw*RQ +(1-specificity)*deltar

        ddeath<-dailydeaths
    return(list(c(ds,de,di,dr,dsq,deq,diq,drq, ddeath)))})}

solveOde<-function(x, parameters){
    # parameters - named vector of parameters
    return(ode(y=State, Mod_times, func=QModel, parms=parameters, method="rk4"))
}


QModelOut<-function(times, X, pars){
    pars_est<- pars

    pars_est["cases_report"] <- 7

    # solve<-ConvertToArray(parLapply(cl=cl, X=1:nrep, fun=solveOde, parameters=pars))
    # solve<-apply(solve, c(1,2), mean)
    solve<-solveOde(parameters = pars)

    solve[,"Death"]<-floor(solve[,"Death"])

    cases_solve_time<-round(solve[,"time"]+pars_est["cases_report"], digits=ndigits)
    cases_solve<-rowSums( solve[which(is.wholenumber(cases_solve_time)),c("SQ","EQ","IQ","RQ")] )

    death_solve_time<-round(solve[,"time"]+death_report, digits=ndigits)#Pars["death_report"
    death_solve<-solve[which(is.wholenumber(death_solve_time)),"Death"]

    cases_solve_rep<-data.frame(rep_time=cases_solve_time[which(is.wholenumber(cases_solve_time))], cases=cases_solve)
    death_solve_rep<-data.frame(rep_time = death_solve_time[which(is.wholenumber(death_solve_time))], death=death_solve )
    
    return(list(solve[which(is.wholenumber(times)),], cases_solve_rep, death_solve_rep))
}

QModelOut2<-function(pars){
    return(QModelOut(Mod_times, State, pars)[[1]])}


QModelCost<-function(Pars){
    # print(Pars, digits=15)
    # solve<-ConvertToArray(parLapply(cl=cl, X=1:nrep, fun=solveOde, parameters=Pars))
    # solve<-apply(solve, c(1,2), mean)
    Pars["cases_report"]<-7

    solve<-solveOde(parameters = Pars)

    solve[,"Death"]<-floor(solve[,"Death"])

    cases_solve_time<-round(solve[,"time"]+Pars["cases_report"], digits=ndigits)
    cases_solve<-rowSums( solve[which(is.wholenumber(cases_solve_time)),c("SQ","EQ","IQ","RQ")] )

    death_solve_time<-round(solve[,"time"]+death_report, digits=ndigits)
    death_solve<-solve[which(is.wholenumber(death_solve_time)),"Death"]
    rm(solve)
    
    mod_cases<-data.frame(time=cases_solve_time[which(is.wholenumber(cases_solve_time))], cases = cases_solve)
    cost_cases<-modCost(model=mod_cases, obs = obs_cases, weight = "std")

    mod_deaths<-data.frame(time=death_solve_time[which(is.wholenumber(death_solve_time))], deaths = death_solve)
    cost_deaths<-modCost(model = mod_deaths, obs = obs_deaths, cost = cost_cases, weight = "std")
    # print(cost_deaths$minlogp)
    # print(cost_deaths$var)

    cost_deaths$residuals %>%
        group_by(name) %>%
        mutate(sd = sd(obs)) -> cost_sd
    # print( sum(log(dnorm(cost_sd$mod, mean = cost_sd$obs, sd = sd(cost_sd$sd)))))
    return(cost_deaths)
}

State<-c(S=33300000,E=0,I=0,R=0,SQ=0,EQ=0,IQ=0,RQ=0,Death = 0)

Mod_times<- seq(0,as.double(as.Date(end_date) - as.Date(start_date)), stepsize)        

death_report <- 1

# read in rds file
mcmc.fit<- readRDS("results/results_20211030/petar/adaptive_MCMC_fit_red_object_nrep_1nbin_0stepsize_0.1_mean_dist_2021-10-15-8.rds")

# run sensitivity

sv<-sensRange(parms = mcmc.fit$bestpar, parInput = mcmc.fit$par, f = QModelOut2, num = 6000)

pdf(paste0("results/adaptive_MCMC_fit_sensitivity.pdf"))
    plot(sv)
    sumsv<-summary(sv)
    plot(sumsv,  xlab="Days since initial infection", ylab="Out of hospital susceptible population", which="S", main="")
    plot(sumsv,  xlab="Days since initial infection", ylab="Out of hospital exposed population", which="E", main="")
    plot(sumsv,  xlab="Days since initial infection", ylab="Out of hospital infected population", which="I", main="")
    plot(sumsv,  xlab="Days since initial infection", ylab="Out of hospital recovered population", which="R", main="")
    plot(sumsv,  xlab="Days since initial infection", ylab="Hospitalised susceptible population", which="SQ", main="")
    plot(sumsv,  xlab="Days since initial infection", ylab="Hospitalised exposed population", which="EQ", main="") 
    plot(sumsv,  xlab="Days since initial infection", ylab="Hospitalised infected population", which="IQ", main="") 
    plot(sumsv,  xlab="Days since initial infection", ylab="Hospitalised recovered population", which="RQ", main="") 
    plot(sumsv,  xlab="Days since initial infection", ylab="COVID-19 deaths", which="Death", main="") 

dev.off()

write.csv(sumsv, paste0("results/sensrange_summary.csv"))