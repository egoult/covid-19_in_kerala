##
## Produces all the variations on the hospitalisation model for comparison
## And returns a csv with output
## Author: Elizabeth Goult
## Date: 07/07/2020
## Comments: Models
##          Original quarantine model
##          Reduced testing rate
##          No travel restrictions
##          No out-of-hospital measures
##          No in hospital quarantining
##          R0=3
##          R0=3 and no out-of-hospital measures


# clear work space
# clear graphics
rm(list=ls())
graphics.off()

# library
require(deSolve)

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
nrep<-30

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

    
NoTravel<-function(time, X, pars){
    # X compartments 
    # S,E,I,R
    # SQ, EQ, IQ, RQ
    # No reduction in travel into the state

    with(as.list(c(X, pars)),{
        H<-sum(c(S,E,I,R))

        w<- 0.000
        omega <-1 
        omegaw<-1
        specificity<- 1
        sensitivity<-0.85
        
        total_delta<-round(rnorm(1, 46000, 2000))
        ninf<-rbinom(1, size=total_delta, prob=global[floor(time)+1]/7.8e+9) 

        deltas<-total_delta-ninf
        deltae<-runif(1,0, ninf)
        deltai<-runif(1,0, ninf  - deltae)
        deltar<-ninf  - deltae - deltai


        if(time<(as.Date("2020-03-23") - as.Date(start_date))){
            lambda<-lambda1
        }else if (time<(as.Date("2020-04-20") - as.Date(start_date))){ 
            lambda<-lambda2
        }else if (time<(as.Date("2020-04-24") - as.Date(start_date))){ 
            lambda<-lambda3
        }else{
            lambda<-lambda4}

        #Non hospitalised compartments
        ds <- -lambda*S*I/H + w*R + omega*SQ + specificity*deltas-total_delta
        de <-  lambda*S*I/H - p*E + (1-sensitivity)*deltae 
        di <-  p*E - r*I - sigma*I + (1-sensitivity)*deltai
        dr <-  r*I + omegaw*RQ - w*R + specificity*deltar

        # hospitalised
        # work out number of deaths
        
        dailydeaths<-rbinom(1,size=round(IQ), prob=d)
        
        dsq <- (1-specificity)*deltas- omega*SQ 
        deq <- sensitivity*deltae- p*EQ 
        diq <- p*EQ + sigma*I + sensitivity*deltai -r*IQ - dailydeaths
        drq <- r*IQ - omegaw*RQ +(1-specificity)*deltar

        ddeath<-dailydeaths
        dnew_cases<-sum(c(p*E + (1-sensitivity)*deltai, p*EQ + sensitivity*deltai))
        dnew_hosp<-sum(c((1-specificity)*deltas, sensitivity*deltae, sigma*I + sensitivity*deltai, (1-specificity)*deltar))
    return(list(c(ds,de,di,dr,dsq,deq,diq,drq, ddeath, dnew_cases, dnew_hosp)))})}


NoQuarantining<-function(time, X, pars){
    # X compartments 
    # SN,EN,IN,RN
    # SP, EP, IP, RP
    # No quarantining of hospitalised population - i.e full population mixing


    with(as.list(c(X, pars)),{
        HN<-sum(c(SN,EN,IN,RN))
        HP<-sum(c(SP,EP,IP,RP))

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

        #Non tested compartments
        dsn <- -lambda*SN*(IN+IP)/(HN+HP) + w*RN + omega*SP + specificity*deltas-total_delta
        den <-  lambda*SN*(IN+IP)/(HN+HP) - p*EN + (1-sensitivity)*deltae 
        din <-  p*EN - r*IN - sigma*IN + (1-sensitivity)*deltai
        drn <-  r*IN + omegaw*RP - w*RN + specificity*deltar

        # tested
        # work out number of deaths
        dailydeaths<-rbinom(1,size=round(IP), prob=d)
        
        dsp <- (1-specificity)*deltas- omega*SP  -lambda*SP*(IN+IP)/(HN+HP) 
        dep <- lambda*SP*(IN+IP)/(HN+HP) + sensitivity*deltae- p*EP 
        dip <- p*EP + sigma*IN + sensitivity*deltai -r*IP - dailydeaths
        drp <- r*IP - omegaw*RP  + (1-specificity)*deltar

        ddeath<-dailydeaths
        dnew_cases<-sum(c(p*EN + (1-sensitivity)*deltai, p*EP + sensitivity*deltai))
        dnew_hosp<-sum(c((1-specificity)*deltas, sensitivity*deltae, sigma*IN + sensitivity*deltai, (1-specificity)*deltar))
    return(list(c(dsn,den,din,drn,dsp,dep,dip,drp, ddeath, dnew_cases, dnew_hosp)))})}

    All<-function(time, X, pars){
    # X compartments 
    # SN,EN,IN,RN
    # SP, EP, IP, RPz
    # No quarantining of hospitalised - i.e full population mixing
    # reduce testing by changing a


    with(as.list(c(X, pars)),{
        HN<-sum(c(SN,EN,IN,RN))
        HP<-sum(c(SP,EP,IP,RP))


        w<- 0.000
        omega <-1
        omegaw<-1
        specificity<- 1
        sensitivity<-0.85
        
        total_delta<-round(rnorm(1, 46000, 2000))
        ninf<-rbinom(1, size=total_delta, prob=global[floor(time)+1]/7.8e+9) 

        deltas<-total_delta-ninf
        deltae<-runif(1,0, ninf)
        deltai<-runif(1,0, ninf  - deltae)
        deltar<-ninf  - deltae - deltai

        if(time<(as.Date("2020-03-23") - as.Date(start_date))){
            lambda<-lambda1
        }else if (time<(as.Date("2020-04-20") - as.Date(start_date))){ 
            lambda<-lambda2
        }else if (time<(as.Date("2020-04-24") - as.Date(start_date))){ 
            lambda<-lambda3
        }else{
            lambda<-lambda4}

        #Non tested compartments
        dsn <- -lambda*SN*(IN+IP)/(HN+HP) + w*RN + omega*SP + specificity*a*deltas+(1-a)*deltas-total_delta
        den <-  lambda*SN*(IN+IP)/(HN+HP) - p*EN + (1-sensitivity)*a*deltae + (1-a)*deltae 
        din <-  p*EN - r*IN - sigma*IN + (1-sensitivity)*a*deltai + (1-a)*deltai
        drn <-  r*IN + omegaw*RP - w*RN + a*specificity*deltar + (1-a)*deltar

        # tested
        # work out number of deaths
        dailydeaths<-rbinom(1,size=round(IP), prob=d)
        
        dsp <- (1-specificity)*a*deltas- omega*SP -lambda*SP*(IN+IP)/(HN+HP) 
        dep <- lambda*SP*(IN+IP)/(HN+HP) + sensitivity*a*deltae- p*EP 
        dip <- p*EP + sigma*IN + sensitivity*a*deltai -r*IP - dailydeaths
        drp <- r*IP - omegaw*RP + a*(1-specificity)*deltar

        ddeath<-dailydeaths
        dnew_cases<-sum(c(p*EN + (1-sensitivity)*a*deltai + (1-a)*deltai, p*EP + sensitivity*a*deltai))
        dnew_hosp<-sum(c((1-specificity)*a*deltas, sensitivity*a*deltae, sigma*IN + sensitivity*a*deltai, (1-specificity)*a*deltar))
    return(list(c(dsn,den,din,drn,dsp,dep,dip,drp, ddeath, dnew_cases, dnew_hosp)))})}

    NRMSD<-function(mod, obs){
        sqdiff<-(mod-obs)^2
        rootmn<-sqrt(mean(sqdiff))
        norm<-rootmn/mean(obs)
        return(norm)}


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

Reporting<-function(solve, hosp, total){
    Active_total_cases<-rowSums(solve[,total])
    Cumulative_total_cases<-solve[,"New"]
    print(paste("Total cases", max(Active_total_cases), max(Cumulative_total_cases)))

    Active_hospital_cases<-rowSums(solve[,hosp])
    Cumulative_hospital_cases<-solve[,"New_hosp"]

    shifted_time<-solve[,"time"]+reporting_delay
    shifted<-which(shifted_time<=max(solve[,"time"]))
    print(paste("Reporting shifted hospitalised cases", max(Active_hospital_cases[shifted]), max(Cumulative_hospital_cases[shifted])))

    Deaths<-solve[,"Death"]
    print(paste0("Deaths ",max(Deaths)))

    print(paste("mortality rates", max(Deaths)/max(Cumulative_hospital_cases[shifted]), max(Deaths)/max(Cumulative_total_cases) ))
}

# end functions

# All model inputs
State<-c(S=33300000,E=0,I=0,R=0,SQ=0,EQ=0,IQ=0,RQ=0,Death = 0, New=0, New_hosp=0)
Mod_times<- seq(0,as.double(as.Date(end_date) - as.Date(start_date)), stepsize) 
parameters<-c(lambda1 = 1.12999975, lambda2=0.10999994, lambda3 = 1.09000037, lambda4 = 1.09000037, w= 0.000, omega = 1, p=0.19999992, r=0.07142859, 
            sigma=0.44220802, omegaw=1, d=0.0004800 , specificity = 1, sensitivity=0.85)


# current model
quarantine_solve<-replicate(n=nrep, ode(y=State,Mod_times,func=QModel, parms=parameters, method="rk4") )
quarantine_solve<-apply(quarantine_solve, c(1,2), mean)#[which(is.wholenumber(times)),]
quarantine_solve_death<-floor(quarantine_solve[,"Death"])[which(is.wholenumber(Mod_times))]

print("reference")
Reporting(quarantine_solve, c("SQ", "EQ", "IQ", "RQ"), c("I","IQ"))


#Reduced testing rate
testing_parameters<-c(parameters,a=0.1)
testing_solve<-replicate(n=nrep, ode(y=State,Mod_times,func=ReducedQModel, parms=testing_parameters, method="rk4") )
testing_solve<-apply(testing_solve, c(1,2), mean)
testing_solve_death<-floor(testing_solve[,"Death"])[which(is.wholenumber(Mod_times))]

print("Reduced testing")
Reporting(testing_solve, c("SQ", "EQ", "IQ", "RQ"), c("I","IQ"))


#No travel cut off 
travel_solve<-replicate(n=nrep, ode(y=State,Mod_times,func=NoTravel, parms=parameters, method="rk4") )
travel_solve<-apply(travel_solve, c(1,2), mean)#[which(is.wholenumber(times)),]
travel_solve_death<-floor(travel_solve[,"Death"])[which(is.wholenumber(Mod_times))]

print("Travel bans")
Reporting(travel_solve, c("SQ", "EQ", "IQ", "RQ"), c("I","IQ"))


# No lockdown
lockdown_parameters<-parameters
lockdown_parameters["lambda2"]<-parameters["lambda1"]#*0.9
lockdown_parameters["lambda3"]<-parameters["lambda1"]#*0.9
lockdown_parameters["lambda4"]<-parameters["lambda1"]#*0.9

lockdown_solve<-replicate(n=nrep, ode(y=State,Mod_times,func=QModel, parms=lockdown_parameters, method="rk4") )
lockdown_solve<-apply(lockdown_solve, c(1,2), mean)#[which(is.wholenumber(times)),]
lockdown_solve_death<-floor(lockdown_solve[,"Death"])[which(is.wholenumber(Mod_times))]
# lockdown_hosp_solve<-lockdown_solve[which(is.wholenumber(Mod_times)),"New_hosp"]

print("No lockdown")
Reporting(lockdown_solve, c("SQ", "EQ", "IQ", "RQ"), c("I","IQ"))


# Shifted lambda
actual_r0<-parameters["lambda1"]/(parameters["r"]+parameters["sigma"])

shift_parameters<-parameters
shift_parameters["lambda1"]<-parameters["lambda1"]*3/actual_r0
shift_parameters["lambda2"]<-parameters["lambda2"]*3/actual_r0
shift_parameters["lambda3"]<-parameters["lambda3"]*3/actual_r0
shift_parameters["lambda4"]<-parameters["lambda4"]*3/actual_r0

shift_solve<-replicate(n=nrep, ode(y=State,Mod_times,func=QModel, parms=shift_parameters, method="rk4") )
shift_solve<-apply(shift_solve, c(1,2), mean)#[which(is.wholenumber(times)),]
shift_solve_death<-floor(shift_solve[,"Death"])[which(is.wholenumber(Mod_times))]

print(" R0 at 3")
Reporting(shift_solve, c("SQ", "EQ", "IQ", "RQ"), c("I","IQ"))


# Shifted lambda, with lockdown removed
shift_lock_parameters<-parameters
shift_lock_parameters["lambda1"]<-parameters["lambda1"]*3/actual_r0
shift_lock_parameters["lambda2"]<-parameters["lambda1"]*3/actual_r0
shift_lock_parameters["lambda3"]<-parameters["lambda1"]*3/actual_r0
shift_lock_parameters["lambda4"]<-parameters["lambda1"]*3/actual_r0

shift_lock_solve<-replicate(n=nrep, ode(y=State,Mod_times,func=QModel, parms=shift_lock_parameters, method="rk4") )
shift_lock_solve<-apply(shift_lock_solve, c(1,2), mean)#[which(is.wholenumber(times)),]
shift_lock_solve_death<-floor(shift_lock_solve[,"Death"])[which(is.wholenumber(Mod_times))]


print("No lockdown, R0 at 3")
Reporting(shift_lock_solve, c("SQ", "EQ", "IQ", "RQ"), c("I","IQ"))



# No quarantine
Noq_state<-c(SN=33300000,EN=0,IN=0,RN=0,SP=0,EP=0,IP=0,RP=0,Death = 0, New=0, New_hosp=0)
Noq_solve<-replicate(n=nrep, ode(y=Noq_state,Mod_times,func=NoQuarantining, parms=parameters, method="rk4") )
Noq_solve<-apply(Noq_solve, c(1,2), mean)
Noq_solve_death<-floor(Noq_solve[,"Death"])[which(is.wholenumber(Mod_times))]



print("No quarantine")
Reporting(Noq_solve, c("SP", "EP", "IP", "RP"), c("IP", "IN"))

# No quarantine, lockdown or travel restrictions

All_state<-c(SN=33300000,EN=0,IN=0,RN=0,SP=0,EP=0,IP=0,RP=0,Death = 0, New=0, New_hosp=0)
All_solve<-replicate(n=nrep, ode(y=All_state, Mod_times,func=All, parms=c(lockdown_parameters, a=0), method="rk4") )
All_solve<-apply(All_solve, c(1,2), mean)#[which(is.wholenumber(times)),]
All_solve_death<-floor(All_solve[,"Death"])[which(is.wholenumber(Mod_times))]


print("All")
Reporting(All_solve, c("SP", "EP", "IP", "RP"), c("IP", "IN"))

# Results
if (!dir.exists('results/')) {
  dir.create('results/')
}

# plot the number of cases
pdf("results/quarantine_cases_compare_plot.pdf")
# ongoing
plot(x=Mod_times, y=log(quarantine_solve[,"I"]+quarantine_solve[,"IQ"]), xlab="Days since initial infection", ylab="Ongoing COVID-19 infections", ylim=c(0, max(log(All_solve[,"IP"]+All_solve[,"IN"]),na.rm=T)), typ="l", yaxt="n", cex.lab=1.2)
axis(2, labels = c(1,10,100,1000,10000,100000, 1000000, 10000000) , at=log(c(1,10,100,1000,10000,100000,1000000,10000000)))
lines(Mod_times, log(testing_solve[,"I"]+testing_solve[,"IQ"]),col=2)
lines(Mod_times, log(travel_solve[,"I"]+travel_solve[,"IQ"]), col=3)
lines(Mod_times, log(lockdown_solve[,"I"]+lockdown_solve[,"IQ"]),  col=4)
lines(Mod_times, log(Noq_solve[,"IN"]+Noq_solve[,"IP"]), col=5)
lines(Mod_times, log(All_solve[,"IN"]+All_solve[,"IP"]), col=6)
lines(Mod_times, log(shift_solve[,"I"]+shift_solve[,"IQ"]),  col=7) # new
lines(Mod_times, log(shift_lock_solve[,"I"]+shift_lock_solve[,"IQ"]),  col=8) # new

legend("topleft", pch=20, col=1:8,  c("Current response","Reduced testing","No travel restrictions","No lock-down","No quarantine", "No response", "R0 of 3", "R0 of 3, no lockdown"), bty="n", cex=1)

log_ongoing_cases<-data.frame(Time = Mod_times, 
                    kerala=log(quarantine_solve[,"I"]+quarantine_solve[,"IQ"]), 
                    testing=log(testing_solve[,"I"]+testing_solve[,"IQ"]),
                    travel=log(travel_solve[,"I"]+travel_solve[,"IQ"]),
                    lockdown = log(lockdown_solve[,"I"]+lockdown_solve[,"IQ"]), 
                    no_quarantine= log(Noq_solve[,"IN"]+Noq_solve[,"IP"]),
                    All = log(All_solve[,"IN"]+All_solve[,"IP"]),
                    shift = log(shift_solve[,"I"]+shift_solve[,"IQ"]), 
                    shifted_lockdown = log(shift_lock_solve[,"I"]+shift_lock_solve[,"IQ"]))
write.csv(log_ongoing_cases, "results/log_ongoing_cases.csv")

ongoing_cases<-data.frame(Time = Mod_times, 
                    kerala=(quarantine_solve[,"I"]+quarantine_solve[,"IQ"]), 
                    testing=(testing_solve[,"I"]+testing_solve[,"IQ"]),
                    travel=(travel_solve[,"I"]+travel_solve[,"IQ"]),
                    lockdown = (lockdown_solve[,"I"]+lockdown_solve[,"IQ"]), 
                    no_quarantine= (Noq_solve[,"IN"]+Noq_solve[,"IP"]),
                    All = (All_solve[,"IN"]+All_solve[,"IP"]),
                    shift = (shift_solve[,"I"]+shift_solve[,"IQ"]), 
                    shifted_lockdown = (shift_lock_solve[,"I"]+shift_lock_solve[,"IQ"]))
write.csv(ongoing_cases, "results/ongoing_cases.csv")



# ongoing vs population size
plot(x=Mod_times, y=(quarantine_solve[,"I"]+quarantine_solve[,"IQ"])/rowSums(quarantine_solve[,2:9]),
     xlab="Days since initial infection", ylab="Ongoing COVID-19 incidence ", 
     ylim=c(0, max((All_solve[,"IP"]+All_solve[,"IN"])/rowSums(All_solve[,2:9]),na.rm=T)), typ="l", cex.lab=1.2)#, yaxt="n")
lines(Mod_times, (testing_solve[,"I"]+testing_solve[,"IQ"])/rowSums(testing_solve[,2:9]),col=2)
lines(Mod_times, (travel_solve[,"I"]+travel_solve[,"IQ"])/rowSums(travel_solve[,2:9]), col=3)
lines(Mod_times, (lockdown_solve[,"I"]+lockdown_solve[,"IQ"])/rowSums(lockdown_solve[,2:9]),  col=4)
lines(Mod_times, (Noq_solve[,"IN"]+Noq_solve[,"IP"])/rowSums(Noq_solve[,2:9]), col=5)
lines(Mod_times, (All_solve[,"IN"]+All_solve[,"IP"])/rowSums(All_solve[,2:9]), col=6)
lines(Mod_times, (shift_solve[,"I"]+shift_solve[,"IQ"])/rowSums(shift_solve[,2:9]),  col=7)#new
lines(Mod_times, (shift_lock_solve[,"I"]+ shift_lock_solve[,"IQ"])/rowSums(shift_lock_solve[,2:9]),  col=8) #new

legend("topleft", pch=20, col=1:8,  c("Current response","Reduced testing","No travel restrictions","No lock-down","No quarantine", "No response", "R0 of 3", "R0 of 3, no lockdown"), bty="n", cex=1)


ongoing_incidence<-data.frame(Time = Mod_times,
                        kerala=(quarantine_solve[,"I"]+quarantine_solve[,"IQ"])/rowSums(quarantine_solve[,2:9]),
                        testing=(testing_solve[,"I"]+testing_solve[,"IQ"])/rowSums(testing_solve[,2:9]),
                        travel=(travel_solve[,"I"]+travel_solve[,"IQ"])/rowSums(travel_solve[,2:9]),
                        lockdown=(lockdown_solve[,"I"]+lockdown_solve[,"IQ"])/rowSums(lockdown_solve[,2:9]),
                        no_quarantine=(Noq_solve[,"IN"]+Noq_solve[,"IP"])/rowSums(Noq_solve[,2:9]),
                        All=(All_solve[,"IN"]+All_solve[,"IP"])/rowSums(All_solve[,2:9]),
                        shift=(shift_solve[,"I"]+shift_solve[,"IQ"])/rowSums(shift_solve[,2:9]),
                        shifted_lockdown=(shift_lock_solve[,"I"]+shift_lock_solve[,"IQ"])/rowSums(shift_lock_solve[,2:9]))

write.csv(ongoing_incidence,"results/ongoing_incidence.csv")

# cumulative
plot(x=Mod_times, y=log(quarantine_solve[,"New"]), xlab="Days since initial infection", ylab="Cumulative COVID-19 infections", ylim=c(0, max(log(All_solve[,"New"]),na.rm=T)), typ="l", yaxt="n", cex.lab=1.2)
axis(2, labels = c(1,10,100,1000,10000,100000, 1000000, 10000000) , at=log(c(1,10,100,1000,10000,100000,1000000,10000000)))
lines(Mod_times, log(testing_solve[,"New"]),col=2)
lines(Mod_times, log(travel_solve[,"New"]), col=3)
lines(Mod_times, log(lockdown_solve[,"New"]),  col=4)
lines(Mod_times, log(Noq_solve[,"New"]), col=5)
lines(Mod_times, log(All_solve[,"New"]), col=6)
lines(Mod_times, log(shift_solve[,"New"]),  col=7) #new
lines(Mod_times, log(shift_lock_solve[,"New"]),  col=8) #new

legend("topleft", pch=20, col=1:8,  c("Current response","Reduced testing","No travel restrictions","No lock-down","No quarantine", "No response", "R0 of 3", "R0 of 3, no lockdown"), bty="n", cex=1)

log_cumulative_cases<-data.frame(Time=Mod_times,
                        kerala=log(quarantine_solve[,"New"]),
                        testing=log(testing_solve[,"New"]),
                        travel=log(travel_solve[,"New"]),
                        lockdown=log(lockdown_solve[,"New"]),
                        no_quarantine=log(Noq_solve[,"New"]),
                        All=log(All_solve[,"New"]),
                        shift=log(shift_solve[,"New"]),
                        shifted_lockdown=log(shift_lock_solve[,"New"]))

write.csv(log_cumulative_cases,"results/log_cumulative_cases.csv")


cumulative_cases<-data.frame(Time=Mod_times,
                        kerala=(quarantine_solve[,"New"]),
                        testing=(testing_solve[,"New"]),
                        travel=(travel_solve[,"New"]),
                        lockdown=(lockdown_solve[,"New"]),
                        no_quarantine=(Noq_solve[,"New"]),
                        All=(All_solve[,"New"]),
                        shift=(shift_solve[,"New"]),
                        shifted_lockdown=(shift_lock_solve[,"New"]))

write.csv(cumulative_cases,"results/cumulative_cases.csv")

# cumulative vs poulation size
plot(x=Mod_times, y=(quarantine_solve[,"New"]/rowSums(quarantine_solve[,2:9])), xlab="Days since initial infection", ylab="Cumulative COVID-19 incidence", ylim=c(0, max((All_solve[,"New"]/rowSums(All_solve[,2:9])),na.rm=T)), typ="l", cex.lab=1.2)#, yaxt="n")
lines(Mod_times, (testing_solve[,"New"]/rowSums(testing_solve[,2:9])),col=2)
lines(Mod_times, (travel_solve[,"New"]/rowSums(travel_solve[,2:9])), col=3)
lines(Mod_times, (lockdown_solve[,"New"]/rowSums(lockdown_solve[,2:9])),  col=4)
lines(Mod_times, (Noq_solve[,"New"]/rowSums(Noq_solve[,2:9])), col=5)
lines(Mod_times, (All_solve[,"New"]/rowSums(All_solve[,2:9])), col=6)
lines(Mod_times, (shift_solve[,"New"]/rowSums(shift_solve[,2:9])),  col=7) #new
lines(Mod_times, (shift_lock_solve[,"New"]/rowSums(shift_lock_solve[,2:9])),  col=8) #new

legend("topleft", pch=20, col=1:8,  c("Current response","Reduced testing","No travel restrictions","No lock-down","No quarantine", "No response", "R0 of 3", "R0 of 3, no lockdown"), bty="n", cex=1)

cumulative_incidence<-data.frame(Time=Mod_times,
                            kerala=(quarantine_solve[,"New"]/rowSums(quarantine_solve[,2:9])),
                            testing=(testing_solve[,"New"]/rowSums(testing_solve[,2:9])),
                            travel=(travel_solve[,"New"]/rowSums(travel_solve[,2:9])),
                            lockdown=(lockdown_solve[,"New"]/rowSums(lockdown_solve[,2:9])),
                            no_quarantine=(Noq_solve[,"New"]/rowSums(Noq_solve[,2:9])),
                            All=(All_solve[,"New"]/rowSums(All_solve[,2:9])),
                            shift=(shift_solve[,"New"]/rowSums(shift_solve[,2:9])),
                            shifted_lockdown=(shift_lock_solve[,"New"]/rowSums(shift_lock_solve[,2:9])) )

write.csv(cumulative_incidence, "results/cumulative_incidence.csv")

dev.off()


pdf("results/supplementary_variations_cases.pdf")
plot(x=Mod_times, y=log(quarantine_solve[,"I"]+quarantine_solve[,"IQ"]), xlab="Days since initial infection", ylab="Ongoing COVID-19 infections", ylim=c(0, max(log(shift_lock_solve[,"I"]+shift_lock_solve[,"IQ"]),na.rm=T)), typ="l", yaxt="n", cex.lab=1.2)
axis(2, labels = c(1,10,100,1000,10000,100000, 1000000, 10000000) , at=log(c(1,10,100,1000,10000,100000,1000000,10000000)))
lines(Mod_times, log(lockdown_solve[,"I"]+lockdown_solve[,"IQ"]),  col=2)
#lines(Mod_times, log(All_solve[,"IN"]+All_solve[,"IP"]), col=6)
lines(Mod_times, log(shift_solve[,"I"]+shift_solve[,"IQ"]),  col=3) # new
lines(Mod_times, log(shift_lock_solve[,"I"]+shift_lock_solve[,"IQ"]),  col=4) # new

legend("topleft", pch=20, col=1:4,  c("Current response","No out-of-hospital measures", expression(R[0]~"="~3), expression(R[0]~"="~3~","~"no"~"lock-down")), bty="n", cex=1)

# cumulative
plot(x=Mod_times, y=log(quarantine_solve[,"New"]), xlab="Days since initial infection", ylab="Cumulative COVID-19 infections", ylim=c(0, max(log(shift_lock_solve[,"New"]),na.rm=T)), typ="l", yaxt="n", cex.lab=1.2)
axis(2, labels = c(1,10,100,1000,10000,100000, 1000000, 10000000) , at=log(c(1,10,100,1000,10000,100000,1000000,10000000)))
lines(Mod_times, log(lockdown_solve[,"New"]),  col=2)
lines(Mod_times, log(shift_solve[,"New"]),  col=3) #new
lines(Mod_times, log(shift_lock_solve[,"New"]),  col=4) #new

legend("topleft", pch=20, col=1:4,  c("Reference model","No out-of-hospital measures", expression(R[0]~"="~3), expression(R[0]~"="~3~","~"no"~"lock-down")), bty="n", cex=1)
dev.off()




# plot number of deaths
pdf("results/quarantine_death_compare_plot.pdf")
days<-Mod_times[which(is.wholenumber(Mod_times))]
plot(x=days, y=log(quarantine_solve_death), xlab="Days since initial infection", ylab="Total COVID-19 deaths", ylim=c(0, max(log(All_solve_death),na.rm=T)), pch=20, col=1, yaxt="n", cex.lab=1.2)
axis(2, labels = c(1,10,100,1000,10000,100000, 1000000, 10000000) , at=log(c(1,10,100,1000,10000,100000,1000000,10000000)))
points(days, log(testing_solve_death),pch=20,col=2)
points(days, log(travel_solve_death),pch=20, col=3)
points(days, log(lockdown_solve_death), pch=20, col=4)
points(days, log(Noq_solve_death),pch=20, col=5)
points(days, log(All_solve_death),pch=20, col=6)
points(days, log(shift_solve_death), pch=20, col=7) #new
points(days, log(shift_lock_solve_death), pch=20, col=8) #new 

legend("topleft", pch=20, col=1:8,  c("Current response","Reduced testing","No travel restrictions","No lock-down","No quarantine", "No response", "R0 of 3", "R0 of 3, no lockdown"), bty="n", cex=1)
dev.off()


pdf("results/supplementary_variations_deaths.pdf")
days<-Mod_times[which(is.wholenumber(Mod_times))]
plot(x=days, y=log(quarantine_solve_death), xlab="Days since initial infection", ylab="Total COVID-19 deaths", ylim=c(0, max(log(shift_lock_solve_death),na.rm=T)), pch=20, col=1, yaxt="n", cex.lab=1.2)
axis(2, labels = c(1,10,100,1000,10000,100000, 1000000, 10000000) , at=log(c(1,10,100,1000,10000,100000,1000000,10000000)))
points(days, log(lockdown_solve_death), pch=20, col=2)
points(days, log(shift_solve_death), pch=20, col=3) #new
points(days, log(shift_lock_solve_death), pch=20, col=4) #new 

legend("topleft", pch=20, col=1:4,  c("Reference model","No out-of-hospital measures", expression(R[0]~"="~3), expression(R[0]~"="~3~","~"no"~"lock-down")), bty="n", cex=1)
dev.off()

log_death_compare<-data.frame(Time=days, 
                    kerala=log(quarantine_solve_death),
                    testing=log(testing_solve_death),
                    travel=log(travel_solve_death),
                    lockdown=log(lockdown_solve_death),
                    no_quarantine=log(Noq_solve_death),
                    All=log(All_solve_death),
                    shift=log(shift_solve_death),
                    shifted_lockdown=log(shift_lock_solve_death) )

write.csv(log_death_compare, "results/log_deaths.csv")


death_compare<-data.frame(Time=days, 
                    kerala=(quarantine_solve_death),
                    testing=(testing_solve_death),
                    travel=(travel_solve_death),
                    lockdown=(lockdown_solve_death),
                    no_quarantine=(Noq_solve_death),
                    All=(All_solve_death),
                    shift=(shift_solve_death),
                    shifted_lockdown=(shift_lock_solve_death) )

write.csv(death_compare, "results/death_compare.csv")



# plot total population size
pdf("results/quarantine_population_compare_plot.pdf")
plot(x=Mod_times, y=rowSums(quarantine_solve[,2:9]), xlab="Days since initial infection", ylab="Kerala population", typ="l", col=1, ylim=c(min(rowSums(All_solve[,2:9])), 1.005*max(rowSums(quarantine_solve[,2:9]))), cex.lab=1.2)#, yaxt="n")
lines(Mod_times, (rowSums(testing_solve[,2:9])),col=2)
lines(Mod_times, (rowSums(travel_solve[,2:9])),col=3)
lines(Mod_times, (rowSums(lockdown_solve[,2:9])), col=4)
lines(Mod_times, (rowSums(Noq_solve[,2:9])), col=5)
lines(Mod_times, (rowSums(All_solve[,2:9])), col=6)
lines(Mod_times, (rowSums(shift_solve[,2:9])), col=7) # new
lines(Mod_times, (rowSums(shift_lock_solve[,2:9])), col=8) # new

legend("topleft", pch=20, col=1:8,  c("Current response","Reduced testing","No travel restrictions","No lock-down","No quarantine", "No response", "R0 of 3", "R0 of 3, no lockdown"), bty="n", cex=1)

population<-data.frame(Time=Mod_times, 
                    testing=(rowSums(testing_solve[,2:9])),
                    travel=(rowSums(travel_solve[,2:9])),
                    lockdown=(rowSums(lockdown_solve[,2:9])),
                    no_quarantine=(rowSums(Noq_solve[,2:9])),
                    All=(rowSums(All_solve[,2:9])),
                    shift=(rowSums(shift_solve[,2:9])),
                    shifted_lockdown=(rowSums(shift_lock_solve[,2:9])) )

write.csv(population, file = "results/population_compare.csv")

current<-33300000
plot(x=Mod_times, y=rowSums(quarantine_solve[,2:9])/current, xlab="Days since initial infection", ylab="Relative population", typ="l", col=1, ylim=c(0.99,1.01), cex.lab=1.2)
lines(Mod_times, (rowSums(testing_solve[,2:9]))/current,col=2)
lines(Mod_times, (rowSums(travel_solve[,2:9]))/current,col=3)
lines(Mod_times, (rowSums(lockdown_solve[,2:9]))/current, col=4)
lines(Mod_times, (rowSums(Noq_solve[,2:9]))/current, col=5)
lines(Mod_times, (rowSums(All_solve[,2:9]))/current, col=6)
lines(Mod_times, (rowSums(shift_solve[,2:9]))/current, col=7) # new
lines(Mod_times, (rowSums(shift_lock_solve[,2:9]))/current, col=8) #new

legend("topleft", pch=20, col=1:8,  c("Reference model","Reduced testing","No travel restrictions","No lock-down","No quarantine", "No response", "R0 of 3", "R0 of 3, no lockdown"), bty="n", cex=1)

relative_population<-population/current

write.csv(relative_population, file = "results/relative_population.csv")
dev.off()


