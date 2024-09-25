#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: Sept. 2024
#' title: Build the survey and action data
#' ---
#' Clear workspace and load packages:
rm(list=ls())
library(data.table)
library(stringr)
library(zoo)
library(ggplot2)
rm(list=ls())

########## Survey data ##########
JMCC <- fread("../../Data/jmcc.csv")
CPSR <- fread("../../Data/cpsr.csv")

setnames(JMCC,
         c("trust_1", "trust_2", "trust_5", "legis_1", "legis_2", "legis_4"),
         c("trustFatah", "trustHamas", "trustNone",  "legisFatah", "legisHamas", "legisNone"))

dat <- merge(JMCC, CPSR, by = c("month", "year"), all.x=T, all.y=T)

dat[, date:=as.yearmon(paste0(year, "-", month), "%Y-%B")]
dat <- dat[order(date)]

dat <- subset(dat, select=c(trustHamas,trustFatah,
                            supportHamas,supportFatah,
                            legisHamas, legisFatah,
                            date))

########## Attack data #############
load("../../Data/gtd.rdata")
gtd <- subset(gtd, iyear > 1992 & country %in% c(97,155))
#fill using acosta, both groups attacked in 12/1993
acosta <- fread("../../Data/acosta1993.csv")


### Pull out just Fatah and Hamas ###
Hamas.names <- c("Hamas (Islamic Resistance Movement)", "Hamas")
Fatah.names <-  c("Al-Fatah",
                 "Palestine Liberation Organization (PLO)")

H1 <- gtd[gname %in% Hamas.names | gname2 %in% Hamas.names | gname3 %in% Hamas.names, list(iyear, imonth, nkill)]
H1[,date:=as.yearmon(paste(iyear, "-", imonth, sep=""))]
H1[,attacks := sum(length(iyear)), by=list(date)]
H1[,fatalities := sum(nkill, na.rm=T), by=list(date)]
H1 <- unique(H1[, `:=`(iyear=NULL, imonth=NULL, nkill=NULL)])
H1 <- H1[date>="Jan 1994"]
H1 <- merge(H1, data.table(date=as.yearmon(seq(as.Date("1994/1/1"), as.Date("2018/12/31"), by="month"))),
            by="date", all=TRUE)
H1[,attacks:=ifelse(is.na(attacks),0, attacks)]
H1[,fatalities:=ifelse(is.na(fatalities),0, fatalities)]

F1 <- gtd[gname %in% Fatah.names | gname2 %in% Fatah.names | gname3 %in% Fatah.names, list(iyear, imonth, nkill)]
F1[,date:=as.yearmon(paste(iyear, "-", imonth, sep=""))]
F1[,attacks := sum(length(iyear)), by=list(date)]
F1[,fatalities := sum(nkill, na.rm=T), by=list(date)]
F1 <- unique(F1[, `:=`(iyear=NULL, imonth=NULL, nkill=NULL)])
F1 <- F1[date>="Jan 1994"]
F1 <- merge(F1, data.table(date=as.yearmon(seq(as.Date("1994/1/1"), as.Date("2018/12/31"), by="month"))),
            by="date", all=TRUE)
F1[,attacks:=ifelse(is.na(attacks),0, attacks)]
F1[,fatalities:=ifelse(is.na(fatalities),0, fatalities)]


names(H1)[2:3] <- c("Hattacks", "Hkills")
names(F1)[2:3] <- c("Fattacks", "Fkills")

# dat <- data.table(dat)
dat <- merge(dat, H1, by="date")
dat <- merge(dat, F1, by="date")


#fill using acosta, both groups attacked in 12/1993
dat[, lag.Hattacks:=shift(Hattacks, fill=nrow(acosta[ORGANIZATION=="Hamas"]))]
dat[, lag.Fattacks:=shift(Fattacks, fill=nrow(acosta[ORGANIZATION=="Fatah"]))]
dat[, Hkills0:=Hkills]
dat[, Fkills0:=Fkills]
dat[, Hkills:=shift(Hkills, fill=sum(acosta[ORGANIZATION=="Hamas" ]$KILLED))]
dat[, Fkills:=shift(Fkills, fill=sum(acosta[ORGANIZATION=="Fatah" ]$KILLED))]



attackHist <- with(dat,data.frame(Attacks=c(Hattacks, Fattacks),
                         Actor=rep(c("Hamas", "Fatah"), each=300)))

figureA1 <- ggplot(attackHist )+
    geom_histogram(aes(x=Attacks), fill="purple", bins=30)+
    facet_wrap(~Actor, ncol=2)+
    theme_bw(16)+
    ylab("Number of Months")+
    xlab("Attacks/month")
ggsave(figureA1, file="../../Output/Figures/figureA1.pdf", height=4, width=5)


cat("Summary stats on Fatah and Hamas attacks/per month\n")
print(summary(dat[,.(Fattacks, Hattacks)]))
save(dat, file="../../Data/actionsSetup.Rdata")
