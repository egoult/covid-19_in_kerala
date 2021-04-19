##
## Takes sensrange_summary.csv and splits it into compartments
## SO,EO,IO,RO,SH,EH,IH,RH
## author: Elizabeth Goult
## date:09/07/2020
## comments: 

# clear workspace
# close plots
rm(list=ls())
graphics.off()

# read in data
sumsv<-read.csv("results/christina_results/sensrange_summary_nrep_100nbin_1000stepsize_0.1.csv")

comp_names<-c("S", "E", "I", "R", "SQ", "EQ", "IQ", "RQ", "Death")
t_length<-121

if(!dir.exists("results/christina_results/model_output_SA/")){
    print("pamnts")
    dir.create("results/christina_results/model_output_SA/")
}
for(i in comp_names){
    # find which

    comp_rows<-which(sumsv[,"X"] %in% paste0(i,1:t_length))
    i_df<-sumsv[comp_rows,]
    write.csv(i_df, paste0("results/christina_results/model_output_SA/",i,"_sensitivity_analysis.csv"))
}

SQ<-read.csv("results/christina_results/model_output_SA/SQ_sensitivity_analysis.csv")
EQ<-read.csv("results/christina_results/model_output_SA/EQ_sensitivity_analysis.csv")
IQ<-read.csv("results/christina_results/model_output_SA/IQ_sensitivity_analysis.csv")
RQ<-read.csv("results/christina_results/model_output_SA/RQ_sensitivity_analysis.csv")

I<-read.csv("results/christina_results/model_output_SA/I_sensitivity_analysis.csv")


cases<-I[,-c(1:3)]+IQ[,-c(1:3)]
cases$Time<-I[,"x"]
write.csv(cases, "results/christina_results/cases_uncertainty.csv")

hosp<-SQ[,-c(1:3)]+EQ[,-c(1:3)]+IQ[,-c(1:3)]+RQ[,-c(1:3)]
hosp$Time<-SQ[,"x"]

write.csv(hosp, "results/christina_results/unshifted_hospitalised_cases_uncertainty.csv")

hosp<-SQ[,-c(1:3)]+EQ[,-c(1:3)]+IQ[,-c(1:3)]+RQ[,-c(1:3)]
hosp$Time<-SQ[,"x"]+5

write.csv(hosp, "results/christina_results/hospitalised_cases_uncertainty.csv")

Deaths<-read.csv("results/christina_results/model_output_SA/Death_sensitivity_analysis.csv")
Deaths<-Deaths[,-c(1:2)]
names(Deaths)[1]<-"Time"

write.csv(Deaths, "results/christina_results/death_uncertainty.csv")
