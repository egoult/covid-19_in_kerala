
require(ggplot2)

rm(list=ls())
graphics.off()

casesiiq<-read.csv("results/christina_results/cases_uncertainty.csv")
hosp<-read.csv("results/christina_results/hospitalised_cases_uncertainty.csv")
kerala<-read.csv("RAW_DATA/kerala_covid19_20200712.csv")[1:(1+ as.Date("2020-05-30")-as.Date("2020-01-30") ),]
kerala$Time<-as.Date(kerala$Date,"%d/%m/%y")-as.Date("2020-01-30")

plt<-ggplot(kerala, aes(x=Time, y=Current_cases))+
    geom_point()+
    geom_line(data=casesiiq, aes(x=Time, y=Mean), color="red")+
    #geom_ribbon(data=casesiiq, aes(group=1, ymin=Min, ymax=Max))+
    geom_line(data=hosp, aes(x=Time, y=Mean), color="blue")+
    #geom_ribbon(data=hosp, aes(x=Time, ymin=Min, ymax=Max, alpha=0.1), color="blue")+
    labs(x="Active, hospitalised COVID-19 cases", y="Days since first observed case")+
    theme_classic()

pdf("figure-2a.pdf")
print(plt)
dev.off()