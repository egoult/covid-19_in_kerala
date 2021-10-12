##
##
##
##
##

# clear workspace
# clear plots
rm(list = ls())
graphics.off()

# packages()
require(tidyverse)
require(ggplot2)
library(deSolve)
library(minpack.lm)
library(FME)

today<-Sys.Date()

# cl<-makeCluster(detectCores()-3)

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
nmcmc<-5000
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

AddMeans<-function(df, mp){
    # Function takes data prame of parameters and sets them all to the same format
    # adds in parameters equal to their converged mean when removed
    # df - data frame of parameters
    # mp - vectir of mean parameter values


    param_names<-names(mp)

    mean_names<-param_names[!(param_names %in% names(df))]
    if(length(mean_names)>0){    
        for(name in mean_names){
            df[,name]<-mp[name]
        }
    }
    return(df)
    }

# data
mcmcfit<-list()
# _2021-08-21
mcmcfit$all<-read.csv("results/MCMCfit_parameters_nrep_1nbin_0stepsize_0.1_mean_dist.csv", row.names = "X")
mcmcfit$constant<-read.csv("results/MCMCfit_parameters_nrep_1nbin_0stepsize_0.1_mean_dist_2021-08-21_constant.csv", row.names = "X")
mcmcfit$constant2<-read.csv("results/MCMCfit_parameters_nrep_1nbin_0stepsize_0.1_mean_dist_2021-08-21_constant_2.csv", row.names = "X")
mcmcfit$sigma<-read.csv("results/MCMCfit_parameters_nrep_1nbin_0stepsize_0.1_mean_dist_2021-08-22_constant.csv", row.names = "X")
mcmcfit$pd<-read.csv("results/MCMCfit_parameters_nrep_1nbin_0stepsize_0.1_mean_dist_2021-08-23_constant_2.csv", row.names = "X")
mcmcfit$pd2<-read.csv("results/MCMCfit_parameters_nrep_1nbin_0stepsize_0.1_mean_dist_2021-08-23_constant_3.csv", row.names = "X")

# mean parameter values
mean_values<-c(lambda1 = 1.21257929, 
            lambda2 = 0.21480796,
            lambda3 = 1.16467335,            
            sigma = 0.520264,
            d = 0.0004897,
            p = 0.2574875,
            r = 0.1053689996,
            cases_report = 6.720643)
death_report<-1

mcmcfit <- lapply(mcmcfit, FUN = AddMeans, mp = mean_values)
mcmc_combined <- mcmcfit %>%
    bind_rows()

mcmc_combined$iteration<-1:nrow(mcmc_combined)
mcmc_combined$cases_report[mcmc_combined$cases_report==7] <- 6.7205

mcmc_combined_long<-mcmc_combined %>%
    pivot_longer(cols = !iteration) %>%
    filter(name != "cases_report")
        

# plot the paramters by iteration
plt<- ggplot(data = mcmc_combined_long, aes(x = iteration, y = value, group = name)) +
    geom_line() +
    facet_wrap(~name, scales = "free") +
    theme_classic()+
    labs(x = "Iteration", y = "Parameter value", group = "Parameter")       


pdf("results/parameter_iterations_MCMC.pdf")    
    print(plt)
dev.off()

# sensitivity analysis ----------
State<-c(S=33300000,E=0,I=0,R=0,SQ=0,EQ=0,IQ=0,RQ=0,Death = 0)

Mod_times<- seq(0,as.double(as.Date(end_date) - as.Date(start_date)), stepsize)     

sv<-sensRange(parms = mean_values, parInput = mcmc_combined[,-ncol(mcmc_combined) ], f = QModelOut2, num = 3000)
plot(sv)
sumsv<-summary(sv)
write.csv(sumsv, paste0("results/sensrange_summary_nrep_",30000,"stepsize_",stepsize,"_mean_dist_combined.csv"))

# confidence interval method 
# read in confidence intervals
mcmc_ci<-list()
mcmc_ci$all<-readRDS("results/MCMCfit_object_nrep_1nbin_0stepsize_0.1_mean_dist.rds")
mcmc_ci$constant<-readRDS("results/MCMCfit_object_nrep_1nbin_0stepsize_0.1_mean_dist_2021-08-21_constant.rds")
mcmc_ci$constant2<-readRDS("results/MCMCfit_object_nrep_1nbin_0stepsize_0.1_mean_dist_2021-08-21_constant_2.rds")
mcmc_ci$sigma<-readRDS("results/MCMCfit_object_nrep_1nbin_0stepsize_0.1_mean_dist_2021-08-22_constant.rds")
mcmc_ci$pd<-readRDS("results/MCMCfit_object_nrep_1nbin_0stepsize_0.1_mean_dist_2021-08-23_constant_2.rds")
mcmc_ci$pd2<-readRDS("results/MCMCfit_object_nrep_1nbin_0stepsize_0.1_mean_dist_2021-08-23_constant_3.rds")

mcmc_cis<-lapply(mcmc_ci, FUN = function(x){return(summary(as.mcmc(x$pars))$quantiles)})

# select the most recent confidence intervals

CompareCIs<-function(cis1, cis2){
    # takes two lists of quantiles for the parameters
    # takes the  values from the second list if contained in both

    cis1_rm<-cis1[!(row.names(cis1) %in% row.names(cis2)),] 
    cis_comb<-rbind(cis1_rm, cis2)
    return(cis_comb)
}

mcmc_ci<-CompareCIs(mcmc_cis[[1]], mcmc_cis[[2]]) %>%
    CompareCIs(., mcmc_cis[[3]]) %>%
    CompareCIs(., mcmc_cis[[4]]) %>%
    CompareCIs(., mcmc_cis[[5]]) %>%
    CompareCIs(., mcmc_cis[[6]])


write.csv(mcmc_ci, "results/mcmc_quantiles.csv")
parameter_box<-mcmc_ci[,c(1, 5)]

# sensitivity analysis 2
sv<-sensRange(parms = mean_values, parRange = parameter_box, f = QModelOut2, num = 3000)
sumsv<-summary(sv)
write.csv(sumsv, paste0("results/sensrange_summary_nrep_",30000,"stepsize_",stepsize,"_mean_dist_cis.csv"))