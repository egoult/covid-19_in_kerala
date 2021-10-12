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
data.path<-paste0("results/sensrange_summary_nrep_",30000,"stepsize_",0.1,"_mean_dist_combined.csv")
sumsv<-read.csv("~/Downloads/All/result_A/sensrange_summary_nrep_1nbin_500stepsize_0.1_mean_dist.csv")

comp_names<-c("S", "E", "I", "R", "SQ", "EQ", "IQ", "RQ", "Death")
t_length<-121

storage<-"results/results_A/"
if(!dir.exists("results/results_A/model_output_SA/")){
    print("creating directory")
    dir.create("results/results_A/")
    dir.create("results/results_A/model_output_SA/")
}

for(i in comp_names){
    # find which

    comp_rows<-which(sumsv[,"X"] %in% paste0(i,1:t_length))
    i_df<-sumsv[comp_rows,]
    write.csv(i_df, paste0(storage,"model_output_SA/",i,"_sensitivity_analysis.csv"))
}

SQ<-read.csv(paste0(storage, "model_output_SA/SQ_sensitivity_analysis.csv"))
EQ<-read.csv(paste0(storage, "model_output_SA/EQ_sensitivity_analysis.csv"))
IQ<-read.csv(paste0(storage, "model_output_SA/IQ_sensitivity_analysis.csv"))
RQ<-read.csv(paste0(storage, "model_output_SA/RQ_sensitivity_analysis.csv"))

I<-read.csv(paste0(storage, "model_output_SA/I_sensitivity_analysis.csv"))


cases<-I[,-c(1:3)]+IQ[,-c(1:3)]
cases$Time<-I[,"x"]
write.csv(cases, paste0(storage,"cases_uncertainty.csv"))

hosp<-SQ[,-c(1:3)]+EQ[,-c(1:3)]+IQ[,-c(1:3)]+RQ[,-c(1:3)]
hosp$Time<-SQ[,"x"]

write.csv(hosp, paste0(storage,"unshifted_hospitalised_cases_uncertainty.csv"))

hosp<-SQ[,-c(1:3)]+EQ[,-c(1:3)]+IQ[,-c(1:3)]+RQ[,-c(1:3)]
hosp$Time<-SQ[,"x"]+7

write.csv(hosp, paste0(storage, "hospitalised_cases_uncertainty.csv"))

Deaths<-read.csv(paste0(storage, "model_output_SA/Death_sensitivity_analysis.csv"))
Deaths<-Deaths[,-c(1:2)]
names(Deaths)[1]<-"Time"

write.csv(Deaths, paste0(storage, "death_uncertainty.csv"))
