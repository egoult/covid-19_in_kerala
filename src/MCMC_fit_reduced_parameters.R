##
## Modelling hospitalisation for COVID-19 transmission
## Author: Elizabeth Goult
## Date: 07/08/2021
## Comments: take mean of distributions rather than average of runs.
## try to get a good MCMC
## set p, r  and cases_delay to realistic parameter values


# clear work space
# clear plots
rm(list=ls())
graphics.off()

#library
library(deSolve)
library(minpack.lm)
library(FME)
require(tidyverse)
# require(parallel)

today<-paste0(Sys.Date())
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
nmcmc<-100

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

NRMSD<-function(mod, obs){
    sqdiff<-(mod-obs)^2
    rootmn<-sqrt(mean(sqdiff))
    norm<-rootmn/mean(obs)
    return(norm)}

# End functions

State<-c(S=33300000,E=0,I=0,R=0,SQ=0,EQ=0,IQ=0,RQ=0,Death = 0)

Mod_times<- seq(0,as.double(as.Date(end_date) - as.Date(start_date)), stepsize)            


lambda1_init<-1.2269
lambda2_init<-0.16575
lambda3_init<-1.24
sigma_init<-0.49715
d_init <- 0.00048402

# p_init<-0.22933
# r_init<-0.10119

# cases_report_init<- 7
death_report<-1

pars_init<-c(lambda1=lambda1_init,lambda2=lambda2_init, lambda3=lambda3_init, 
            sigma=sigma_init, d = d_init)#, #p=p_init, r=r_init,
            #cases_report = cases_report_init)

# read in kerala data
kerala<-read.csv("RAW_DATA/kerala_covid19_20200712.csv")
obs_date<-as.Date(kerala$Date, "%d/%m/%y")
cases<-kerala$Current_cases[which(obs_date==start_date):which(obs_date==end_date)] 
deaths<-kerala$Cumulative_deaths[which(obs_date==start_date):which(obs_date==end_date)]

# initial parameters results
# pdf("results/delay_fit_start_mean_dist.pdf")
# start_est<-QModelOut(times=Mod_times, X=State, pars_init)
# start_cases<-start_est[[2]]
# start_death<-start_est[[3]]
#  # start compare
#  print("start")
#  print("cases")
#  print(NRMSD(start_cases[which(start_cases[,"rep_time"] %in% 1:length(cases)),"cases"], cases[which(1:length(cases) %in% start_cases[,"rep_time"])]))
#  print("deaths")
#  print(NRMSD(start_death[which(start_death[,"rep_time"] %in% 1:length(deaths)),"death"], deaths[which(1:length(deaths) %in% start_death[,"rep_time"])]))

# par(mfrow=c(2,1), mar=c(2.5, 4.1, 4.1, 2.1))
#     plot(1:length(cases), cases, main="Initial estimate",ylab="COVID-19 cases", xlim=c(0,length(cases)),ylim=c(0,max(c(cases, start_cases[,"cases"]))) )
#     points(start_cases[1:118,"rep_time"], start_cases[1:118,"cases"], pch=2, col="blue")
#     abline(v=(as.Date("2020-03-23") - as.Date(start_date)), col="red")
#     abline(v=(as.Date("2020-04-20") - as.Date(start_date)), col="red")
#     abline(v=(as.Date("2020-04-24") - as.Date(start_date)), col="red")
#     legend("topright", col=c("black","blue"), pch=c(1,2), c("Observation","Model"), bty="n")

#     par(mar=c(5.1, 4.1,2,2.1))
#     plot(1:length(deaths), deaths, xlab="Days since initial infection", ylab="COVID-19 deaths", xlim=c(0, length(deaths)), ylim=c(0, max(c(deaths, start_death[,"death"]))) )
#     points(start_death[,"rep_time"], start_death[,"death"], pch=2, col="blue")
#     abline(v=(as.Date("2020-03-23") - as.Date(start_date)), col="red")
#     abline(v=(as.Date("2020-04-20") - as.Date(start_date)), col="red")
#     abline(v=(as.Date("2020-04-24") - as.Date(start_date)), col="red")
#     legend("topright", col=c("black","blue"), pch=c(1,2), c("Observation","Model"), bty="n")

# dev.off()

# stop()

# FME fit
cases_weights<-rep(1,length(cases))
obs_cases<-data.frame(time=1:length(cases), cases=cases)
obs_deaths<-data.frame(time = 1:length(deaths), deaths = deaths)


#  LMFit<-modFit(f=QModelCost, p=pars_init, method="Pseudo",lower=c(lambda1=0, lambda2=0, lambda3=0, sigma=0, d=0, p=0, r=0, cases_report=0),
#                upper=c(lambda1=2, lambda2=2, lambda3=2, sigma=1, d=1, p=1, r=2, cases_report=20),
#                control=list(numiter=100)) #500
#  print("LMFit")
#  print(summary(LMFit))

# # save results
# if (!dir.exists('results/')) {
#   dir.create('results/')}

# pdf("results/delay_fit_LM_mean_dist.pdf")
# LM_est<-QModelOut(times=Mod_times, X=State, LMFit$par)
# LM_cases<-LM_est[[2]]
# LM_death<-LM_est[[3]]
 
# par(mfrow=c(2,1), mar=c(2.5, 4.1, 4.1, 2.1))
#     plot(1:length(cases), cases, main="LM estimate",ylab="COVID-19 cases", ylim=c(0,max(c(cases, LM_cases[,"cases"]))) )
#     points(LM_cases[1:(length(cases)-round(LMFit$par["cases_report"])),"rep_time"], LM_cases[1:(length(cases)-round(LMFit$par["cases_report"])),"cases"], pch=2, col="blue")
#     abline(v=(as.Date("2020-03-23") - as.Date(start_date)), col="red")
#     abline(v=(as.Date("2020-04-20") - as.Date(start_date)), col="red")
#     abline(v=(as.Date("2020-04-24") - as.Date(start_date)), col="red")
#     legend("topright", col=c("black","blue"), pch=c(1,2), c("Observation","Model"), bty="n")

#     par(mar=c(5.1, 4.1,2,2.1))
#     plot(1:length(deaths), deaths, xlab="Days since initial infection", ylab="COVID-19 deaths", ylim=c(0, max(c(deaths, LM_death[,"death"]))) )
#     points(LM_death[,"rep_time"], LM_death[,"death"], pch=2, col="blue")
#     abline(v=(as.Date("2020-03-23") - as.Date(start_date)), col="red")
#     abline(v=(as.Date("2020-04-20") - as.Date(start_date)), col="red")
#     abline(v=(as.Date("2020-04-24") - as.Date(start_date)), col="red")
#     legend("topright", col=c("black","blue"), pch=c(1,2), c("Observation","Model"), bty="n")

# dev.off()



# MCMC Fit
# var0<-LMFit$var_ms_unweighted
#cov0<-summary(LMFit)$cov.scaled*0.01
#stop()
tick<- Sys.time()
 MCMCFit<-modMCMC(f =  QModelCost,
                  p = pars_init,
                  jump = 2e-2 * pars_init, 
                  var0 = NULL,
                  wvar0 = NULL,
                  niter=nmcmc,#00, 
                  burninlength=nbin,#00,
                  lower = rep(0, length(pars_init)),
                  upper= c(Inf,Inf,Inf,1,1),#,Inf,Inf,7),
                  verbose=T)
tock<- Sys.time()


print(tock - tick)

 print(summary(MCMCFit))
 print(MCMCFit$bestpar)

 print( summary(as.mcmc(MCMCFit$pars)) )

 plot(MCMCFit)

# plot MCMC fit
pdf(paste0("results/MCMC_fit_reduced_",today,".pdf"))
 mcmc_pars<-c(unlist(c(summary(MCMCFit)["mean",1:length(pars_init)])), cases_report = 7)
 mcmc_est<-QModelOut(times=Mod_times, X=State, mcmc_pars)
 mcmc_cases<-mcmc_est[[2]]
 mcmc_death<-mcmc_est[[3]]
 
 # MCMC compare
 print("MCMC mean")
 print("cases")
 print(NRMSD(mcmc_cases[which(mcmc_cases[,"rep_time"] %in% 1:length(cases)),"cases"], cases[which(1:length(cases) %in% mcmc_cases[,"rep_time"])]))
 print("deaths")
 print(NRMSD(mcmc_death[which(mcmc_death[,"rep_time"] %in% 1:length(deaths)),"death"], deaths[which(1:length(deaths) %in% mcmc_death[,"rep_time"])]))

 best_pars<- c(MCMCFit$bestpar, cases_report = 7)
 best_est<- QModelOut(times=Mod_times, X=State, best_pars)
 best_cases<-best_est[[2]]
 best_death<-best_est[[3]]
 
 # MCMC compare
 print("MCMC best")
 print("cases")
 print(NRMSD(best_cases[which(best_cases[,"rep_time"] %in% 1:length(cases)),"cases"], cases[which(1:length(cases) %in% best_cases[,"rep_time"])]))
 print("deaths")
 print(NRMSD(best_death[which(best_death[,"rep_time"] %in% 1:length(deaths)),"death"], deaths[which(1:length(deaths) %in% best_death[,"rep_time"])]))


 par(mfrow=c(2,1), mar=c(2.5, 4.1, 4.1, 2.1))
     plot(1:length(cases), cases, main="MCMC estimate",ylab="COVID-19 cases", ylim=c(0, max(c(cases,mcmc_cases[,"cases"]))) )
     points(mcmc_cases[1:(length(cases)-round(mcmc_pars["cases_report"])),"rep_time"], mcmc_cases[1:(length(cases)-round(mcmc_pars["cases_report"])),"cases"], pch=2, col="blue")
     abline(v=(as.Date("2020-03-23") - as.Date(start_date)), col="red")
     abline(v=(as.Date("2020-04-20") - as.Date(start_date)), col="red")
     legend("topright", col=c("black","blue"), pch=c(1,2), c("Observation","Model"), bty="n")

     par(mar=c(5.1, 4.1,2,2.1))
     plot(1:length(deaths), deaths, xlab="Days since initial infection", ylab="COVID-19 deaths", ylim=c(0, max(c(deaths, mcmc_death[,"death"]))) )
     points(mcmc_death[,"rep_time"], mcmc_death[,"death"], pch=2, col="blue")
     abline(v=(as.Date("2020-03-23") - as.Date(start_date)), col="red")
     abline(v=(as.Date("2020-04-20") - as.Date(start_date)), col="red")
     legend("topright", col=c("black","blue"), pch=c(1,2), c("Observation","Model"), bty="n")

par(mfrow=c(2,1), mar=c(2.5, 4.1, 4.1, 2.1))

    plot(1:length(cases), cases, main="MCMC best estimate",ylab="COVID-19 cases", ylim=c(0, max(c(cases, best_cases[,"cases"]))) )
     points( best_cases[1:(length(cases)-round(best_pars["cases_report"])),"rep_time"], best_cases[1:(length(cases)-round(best_pars["cases_report"])),"cases"], pch=2, col="blue")
     abline(v=(as.Date("2020-03-23") - as.Date(start_date)), col="red")
     abline(v=(as.Date("2020-04-20") - as.Date(start_date)), col="red")
     legend("topright", col=c("black","blue"), pch=c(1,2), c("Observation","Model"), bty="n")

     par(mar=c(5.1, 4.1,2,2.1))
     plot(1:length(deaths), deaths, xlab="Days since initial infection", ylab="COVID-19 deaths", ylim=c(0, max(c(deaths, best_death[,"death"]))) )
     points(best_death[,"rep_time"], best_death[,"death"], pch=2, col="blue")
     abline(v=(as.Date("2020-03-23") - as.Date(start_date)), col="red")
     abline(v=(as.Date("2020-04-20") - as.Date(start_date)), col="red")
     legend("topright", col=c("black","blue"), pch=c(1,2), c("Observation","Model"), bty="n")

dev.off()


# plot histograms of parameter variables
pdf(paste0("results/MCMC_fit_reduced_histograms_",today,".pdf"))
    hist(MCMCFit, Full=T)

    plot(MCMCFit, Full=T)

    # sv<-sensRange(parms = pars_init, parInput = MCMCFit$par, f = QModelOut2, num = MCMCFit$naccapted)
    # plot(sv)
    # sumsv<-summary(sv)
    # plot(sumsv,  xlab="Days since initial infection", ylab="Out of hospital susceptible population", which="S", main="")
    # plot(sumsv,  xlab="Days since initial infection", ylab="Out of hospital exposed population", which="E", main="")
    # plot(sumsv,  xlab="Days since initial infection", ylab="Out of hospital infected population", which="I", main="")
    # plot(sumsv,  xlab="Days since initial infection", ylab="Out of hospital recovered population", which="R", main="")
    # plot(sumsv,  xlab="Days since initial infection", ylab="Hospitalised susceptible population", which="SQ", main="")
    # plot(sumsv,  xlab="Days since initial infection", ylab="Hospitalised exposed population", which="EQ", main="") 
    # plot(sumsv,  xlab="Days since initial infection", ylab="Hospitalised infected population", which="IQ", main="") 
    # plot(sumsv,  xlab="Days since initial infection", ylab="Hospitalised recovered population", which="RQ", main="") 
    # plot(sumsv,  xlab="Days since initial infection", ylab="COVID-19 deaths", which="Death", main="") 
    


dev.off()
# stopCluster(cl)

write.csv(MCMCFit$par,paste0("results/MCMCfit_reduced_parameters_nrep_",nrep,"nbin_", nbin,"stepsize_",stepsize,"_mean_dist_",today,".csv"))
# write.csv(sumsv, paste0("results/sensrange_summary_nrep_",nrep,"nbin_", nbin,"stepsize_",stepsize,"_mean_dist_",today,".csv"))
saveRDS(MCMCFit, file = paste0("results/MCMCfit_reduced_object_nrep_",nrep,"nbin_", nbin,"stepsize_",stepsize,"_mean_dist_",today,".rds"))

# Regression
# pdf("results/delay_regression_mean_dist.pdf")
#     plot(cases[which(1:length(cases) %in% mcmc_cases[,"rep_time"])], 
#          mcmc_cases[which(mcmc_cases[,"rep_time"] %in% 1:length(cases)),"cases"], 
#          xlab="Observations", ylab="Model", main="Ongoing recorded COVID-19 cases")
#     abline(0,1, col="red")
#     print("cases regression")
#     print(summary(lm(mcmc_cases[which(mcmc_cases[,"rep_time"] %in% 1:length(cases)),"cases"]~cases[which(1:length(cases) %in% mcmc_cases[,"rep_time"])])))

#     plot(deaths[which(1:length(deaths) %in% mcmc_death[,"rep_time"])], 
#          mcmc_death[which(mcmc_death[,"rep_time"] %in% 1:length(deaths)),"death"],
#          xlab="Observations", ylab="Model", main="Total COVID-19 deaths")
#     abline(0,1, col="red")
#     print("deaths regression")
#     print(summary(lm(mcmc_death[which(mcmc_death[,"rep_time"] %in% 1:length(deaths)),"death"]~deaths[which(1:length(deaths) %in% mcmc_death[,"rep_time"])])))
# dev.off()

