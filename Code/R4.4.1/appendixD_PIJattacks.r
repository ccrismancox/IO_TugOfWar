#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: March 2022
#' title: Build the data on PIJ attacks 
#' ---
library(data.table)
library(stringr)
library(zoo)
library(ggplot2)
rm(list=ls())


###############################################
# Attack data
load("../../Data/gtd.rdata")

gtd <- subset(gtd, iyear > 1992 & country %in% c(97,155))


PIJ.names <- "Palestinian Islamic Jihad (PIJ)"
P1 <- gtd[gname %in% PIJ.names | gname2 %in% PIJ.names | gname3 %in% PIJ.names, list(iyear, imonth, nkill)]
P1[,date:=as.yearmon(paste(iyear, "-", imonth, sep=""))]
P1[,attacks := sum(length(iyear)), by=list(date)]
P1[,fatalities := sum(nkill, na.rm=T), by=list(date)]
P1 <- unique(P1[, `:=`(iyear=NULL, imonth=NULL, nkill=NULL)])
P1 <- P1[date>="Jan 1994"]
P1 <- merge(P1, data.table(date=as.yearmon(seq(as.Date("1994/1/1"), as.Date("2018/12/31"), by="month"))),
            by="date", all=TRUE)
P1[,attacks:=ifelse(is.na(attacks),0, attacks)]
P1[,fatalities:=ifelse(is.na(fatalities),0, fatalities)]

names(P1)[2:3] <- c("Pattacks", "Pkills")
other.attack <- P1
save(other.attack, file="../../Data/otherattacks.rdata")
