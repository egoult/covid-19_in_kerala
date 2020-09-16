##
## Takes sensrange_summary.csv and splits it into compartments
## SO,EO,IO,RO,SH,EH,IH,RH
## author:elg@pml.ac.uk
## date:09/07/2020
## comments: fml

# clear worrkspace
# close plots
rm(list=ls())
graphics.off()

# read in data
sumsv<-read.csv("sensrange_summary.csv")

comp_names<-c("S", "E", "I", "R", "SQ", "EQ", "IQ", "RQ", "Death")
t_length<-121

for(i in comp_names){
    # find which

    comp_rows<-which(sumsv[,"X"] %in% paste0(i,1:t_length))
    i_df<-sumsv[comp_rows,]
    write.csv(i_df, paste0("code/model_output_SA/",i,"_sensitivity_analysis.csv"))
}

SQ<-read.csv("SQ_sensitivity_analysis.csv")
EQ<-read.csv("EQ_sensitivity_analysis.csv")
IQ<-read.csv("IQ_sensitivity_analysis.csv")
RQ<-read.csv("RQ_sensitivity_analysis.csv")

hosp<-SQ[,-c(1:3)]+EQ[,-c(1:3)]+IQ[,-c(1:3)]+RQ[,-c(1:3)]
hosp$Time<-SQ[,"x"]+5

write.csv(hosp, "hospitalised_cases_uncertainty.csv")

Deaths<-read.csv("Death_sensitivity_analysis.csv")
Deaths<-Deaths[,-c(1:2)]
names(Deaths)[1]<-"Time"

write.csv(Deaths, "death_uncertainty.csv")
