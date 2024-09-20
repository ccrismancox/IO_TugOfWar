#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: Sept 2024
#' title: Build the survey and action by attack types
#' ---
#' Clear workspace and load packages:
rm(list=ls())
library(data.table)
library(stringr)
library(zoo)
library(ggplot2)

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

###############################################
# Attack data
gtd <- fread("../../Data/gtd.csv")
gtd <- subset(gtd, iyear > 1992 & country %in% c(97,155))
acosta <- fread("../../Data/acosta1993.csv")

Hamas.names <- c("Hamas (Islamic Resistance Movement)", "Hamas")
Fatah.names <-  c("Al-Fatah",
                 "Palestine Liberation Organization (PLO)")

# First pass at separating the attack target codes in a coherent manner
civ.attack.codes <- c(1,5,6,8,9,10, 11, 13,14,15,16,18,19,21)
gov.attack.codes <- c(2,7)
mil.attack.codes <- c(3,4)
nonstate.mil.attack <- c(12,17, 22)
unknown <- c(20)
soft <- c(5,8,10,14,15,18)
soft.sub <- c(2:9, 11:12, 112, 99:102) 
H1 <- gtd[gname %in% Hamas.names | gname2 %in% Hamas.names | gname3 %in% Hamas.names, 
          list(iyear, imonth, nkill, targtype1, targsubtype1)]
H1[,date:=as.yearmon(paste(iyear, "-", imonth, sep=""))]

## This is new ## 
H1[,mil.attacks := sum(targtype1 %in% c(gov.attack.codes,mil.attack.codes)), by=list(date)]
H1[,civ.attacks := sum(targtype1 %in% civ.attack.codes), by=list(date)]
H1[,not.civ.attacks := sum(!targtype1 %in% civ.attack.codes), by=list(date)]
H1[,mil.kill := sum(nkill*(targtype1 %in% c(gov.attack.codes,mil.attack.codes)), na.rm=T), by=list(date)]
H1[,civ.kill  := sum(nkill*(targtype1 %in% civ.attack.codes), na.rm=T), by=list(date)]
H1[,not.civ.kill  := sum(nkill*(!targtype1 %in% civ.attack.codes), na.rm=T), by=list(date)]
H1[,soft:=sum(targtype1 %in% soft | targsubtype1 %in% soft.sub, na.rm=T), by=list(date)]
H1[,targtype1:=NULL]
H1[,targsubtype1:=NULL]
###################

H1[,attacks := sum(length(iyear)), by=list(date)]
H1[,fatalities := sum(nkill, na.rm=T), by=list(date)]

H1 <- unique(H1[, `:=`(iyear=NULL, imonth=NULL, nkill=NULL)])
H1 <- H1[date>="Jan 1994"]
H1 <- merge(H1, data.table(date=as.yearmon(seq(as.Date("1994/1/1"), as.Date("2018/12/31"), by="month"))),
            by="date", all=TRUE)
H1[,attacks:=ifelse(is.na(attacks),0, attacks)]
H1[,fatalities:=ifelse(is.na(fatalities),0, fatalities)]

## This is new ## 
H1[,mil.attacks:=ifelse(is.na(mil.attacks),0, mil.attacks)]
H1[,civ.attacks:=ifelse(is.na(civ.attacks),0, civ.attacks)]
H1[,not.civ.attacks:=ifelse(is.na(not.civ.attacks),0, not.civ.attacks)]
H1[,mil.kill:=ifelse(is.na(mil.kill),0, mil.kill)]
H1[,civ.kill:=ifelse(is.na(civ.kill),0, civ.kill)]
H1[,not.civ.kill:=ifelse(is.na(not.civ.kill),0, not.civ.kill)]
H1[,soft:=ifelse(is.na(soft),0, soft)]

####################



F1 <- gtd[gname %in% Fatah.names | gname2 %in% Fatah.names | gname3 %in% Fatah.names, 
          list(iyear, imonth, nkill, targtype1,targsubtype1)]
F1[,date:=as.yearmon(paste(iyear, "-", imonth, sep=""))]
## This is new ## 
F1[,mil.attacks := sum(targtype1 %in% c(gov.attack.codes,mil.attack.codes)), by=list(date)]
F1[,civ.attacks := sum(targtype1 %in% civ.attack.codes), by=list(date)]
F1[,not.civ.attacks := sum(!targtype1 %in% civ.attack.codes), by=list(date)]
F1[,mil.kill := sum(nkill*(targtype1 %in% c(gov.attack.codes,mil.attack.codes)), na.rm=T), by=list(date)]
F1[,civ.kill  := sum(nkill*(targtype1 %in% civ.attack.codes), na.rm=T), by=list(date)]
F1[,not.civ.kill  := sum(nkill*(!targtype1 %in% civ.attack.codes), na.rm=T), by=list(date)]
F1[,soft:=sum(targtype1 %in% soft | targsubtype1 %in% soft.sub, na.rm=T), by=list(date)]
F1[,targtype1:=NULL]
F1[,targsubtype1:=NULL]

###################
F1[,attacks := sum(length(iyear)), by=list(date)]
F1[,fatalities := sum(nkill, na.rm=T), by=list(date)]
F1 <- unique(F1[, `:=`(iyear=NULL, imonth=NULL, nkill=NULL)])
F1 <- F1[date>="Jan 1994"]
F1 <- merge(F1, data.table(date=as.yearmon(seq(as.Date("1994/1/1"), as.Date("2018/12/31"), by="month"))),
            by="date", all=TRUE)
F1[,attacks:=ifelse(is.na(attacks),0, attacks)]
F1[,fatalities:=ifelse(is.na(fatalities),0, fatalities)]

## This is new ## 
F1[,mil.attacks:=ifelse(is.na(mil.attacks),0, mil.attacks)]
F1[,civ.attacks:=ifelse(is.na(civ.attacks),0, civ.attacks)]
F1[,not.civ.attacks:=ifelse(is.na(not.civ.attacks),0, not.civ.attacks)]
F1[,mil.kill:=ifelse(is.na(mil.kill),0, mil.kill)]
F1[,civ.kill:=ifelse(is.na(civ.kill),0, civ.kill)]
F1[,not.civ.kill:=ifelse(is.na(not.civ.kill),0, not.civ.kill)]
F1[,soft:=ifelse(is.na(soft),0, soft)]

####################

names(H1)[-1] <- c("Hattacks.mil", "Hattacks.civ", "Hattacks.notciv",
                   "Hkills.mil", "Hkills.civ", "Hkills.notciv",
                   "Hsoft",
                   "Hattacks", "Hkills")
names(F1)[-1] <- c("Fattacks.mil", "Fattacks.civ", "Fattacks.notciv", 
                   "Fkills.mil", "Fkills.civ", "Fkills.notciv",
                   "Fsoft",
                   "Fattacks", "Fkills")

# dat <- data.table(dat)
dat <- merge(dat, H1, by="date")
dat <- merge(dat, F1, by="date")




#fill using acosta, both groups attacked in 12/1993
dat[, lag.Hattacks:=shift(Hattacks, fill=nrow(acosta[ORGANIZATION=="Hamas"]))]
dat[, lag.Fattacks:=shift(Fattacks, fill=nrow(acosta[ORGANIZATION=="Fatah"]))]
dat[, lag.Hattacks.mil:=shift(Hattacks.mil, fill=nrow(acosta[ORGANIZATION=="Hamas" & `TARGET TYPE` %in% c("Military", "Police")]))]
dat[, lag.Fattacks.mil:=shift(Fattacks.mil, fill=nrow(acosta[ORGANIZATION=="Fatah" & `TARGET TYPE` =="Military"]))]
dat[, lag.Hattacks.civ:=shift(Hattacks.civ, fill=nrow(acosta[ORGANIZATION=="Hamas" & `TARGET TYPE` =="Private citizens and property"]))]
dat[, lag.Fattacks.civ:=shift(Fattacks.civ, fill=nrow(acosta[ORGANIZATION=="Hamas" & `TARGET TYPE` =="Maritime"]))]
dat[, lag.Hattacks.notciv:=shift(Hattacks.notciv, fill=nrow(acosta[ORGANIZATION=="Hamas" & `TARGET TYPE` !="Private citizens and property"]))]
dat[, lag.Fattacks.notciv:=shift(Fattacks.notciv, fill=nrow(acosta[ORGANIZATION=="Hamas" & `TARGET TYPE` !="Maritime"]))]
dat[, Hkills:=shift(Hkills, fill=sum(acosta[ORGANIZATION=="Hamas" ]$KILLED))]
dat[, Fkills:=shift(Fkills, fill=sum(acosta[ORGANIZATION=="Fatah" ]$KILLED))]
dat[, Hsoft:=shift(Hsoft, fill=1)]
dat[, Fsoft:=shift(Fsoft, fill=0)]

dat.Civs <- dat

save(dat.Civs, file="../../Data/actionsSetup_byAttackType.Rdata")
