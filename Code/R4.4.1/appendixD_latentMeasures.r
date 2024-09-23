#' ---
#' author: Casey Crisman-Cox and Mike Gibilisco
#' date: Sept. 2024
#' title: Measurement models for latent attitudes toward violence and umeployment
#' ---
#' 
#' 
#' Clear workspace and load packages:
rm(list=ls())

######################################
# stuff i need
library(matrixStats)
library(MARSS)
library(data.table)
library(knitr)
library(zoo)
###############################################
# data
JMCC <- fread("../../Data/jmcc.csv")
CPSR <- fread("../../Data/cpsr.csv")

setnames(JMCC, 
         c("trust_1", "trust_2", "trust_5", "legis_1", "legis_2", "legis_4"), 
         c("trustFatah", "trustHamas", "trustNone",  "legisFatah", "legisHamas", "legisNone"))

dat <- merge(JMCC, CPSR, by = c("month", "year"), all.x=T, all.y=T)

dat[, date:=as.yearmon(paste0(year, "-", month), "%Y-%B")]
dat <- dat[order(date)]
###############################################
# recode variables that have missing extremes (strongly support/strongly oppose)

# helper function to accomidate NAs 
`%+na%` <- function(x,y) {ifelse( is.na(x), y, ifelse( is.na(y), x, x+y) )}

dat[,`:=`(pp_sup = supp_1 %+na% supp_2, # support peace process
          pp_opp = supp_3 %+na% supp_4, #oppose peace process 
          bomb_sup = bomb_1 %+na% bomb_2, # support suicide bombings
          bomb_opp =  bomb_3 %+na% bomb_4,#oppose suicide bombings
          intif_sup = intif_1 %+na% intif_2, # support second intifada
          intif_opp = intif_3 %+na% intif_4)] #oppose second intifada


###############################################
# create data sets used to estimate dynamic factor model

dat.employ <- subset(dat, select=c(unemp, unemp_respondent,
                                   unemp_cpsr, 
                                   unemp_pcbs))

dat.violence <- subset(dat, select=c(bination_1, bination_2, bination_3, bination_4,
                                     pp_sup, pp_opp,
                                     milit_1, milit_2,
                                     bomb_sup, bomb_opp,
                                     optp_1,optp_2, optp_3, optp_4,
                                     oslo_1,oslo_2,oslo_3,
                                     intif_sup, intif_opp,
                                     deadp_1, deadp_2, deadp_3,
                                     armedattacks, 
                                     armedattacks_civilians, 
                                     armedattacks_soldiers, armedattacks_settlers))

cntl.listA = list(maxit=2500, minit=250, abstol = 1e-5, conv.test.slope.tol = 0.1)
cntl.listB = list(maxit=25000, minit=1000, abstol = 1e-5, conv.test.slope.tol = 0.1)



#### In text remarks
cat("Minimum observed support for violence generally is",
    min(dat.violence$armedattacks, na.rm=TRUE),
    "percent\n")

cat("Median and IQRs for the different attack targets \n")
round(dat.violence[, lapply(.SD, quantile, c(0.5, .25, .75), na.rm=TRUE),
             .SDcols = c("armedattacks_soldiers", 
                         "armedattacks_settlers",
                         "armedattacks_civilians")]
)





####################################################
# unemployment

dyn_mod_uemp <- list(Z=matrix(names(dat.employ)),
                     A="zero", 
                     R="identity", 
                     B="unconstrained",
                     U="zero",
                     Q="identity", 
                     x0="zero",
                     V0=matrix(3,1,1),
                     tinitx=1)
cntl.list = list(maxit=10000, minit=2000, abstol = 1e-5, conv.test.slope.tol = 0.05)
dyn_mars_uemp <- MARSS(t(scale(dat.employ)), 
                       model=dyn_mod_uemp,
                       control=cntl.list, silent=T)

####################################################
# violence 

dyn_mod_violence1D <- list(Z=matrix(names(dat.violence)),
                           A="zero", 
                           R="identity", 
                           B="unconstrained", 
                           U="zero",
                           Q="identity",
                           x0="zero",
                           V0=matrix(3,1,1),
                           tinitx=1)

dyn_mars_violence1D <- MARSS(t(scale(dat.violence)), 
                             model=dyn_mod_violence1D,
                             control=cntl.list, silent=T)



###############Tables#######################################
marsOut <- data.frame(V1=dyn_mars_uemp$par$Z)
marsOut$Variable <- c("Self reported unemployment rate (JMCC)",
                      "Self reported unemployment rate (PCPSR)",
                      "Estimated unemployment rate (PCPSR)",
                      "Estimated unemployment rate (PCBS)")


marsOut <- marsOut[  ,c("Variable","V1")]
cat(kable(marsOut, row.names=FALSE, col.names=c("Variable", "Est."),
          caption="ML estimates for latent unemployment conditions",
          digits=2, format="pipe"),
    file="../../Output/Tables/tableD2.txt",
    sep="\n")


marsOut <- data.frame(V1=dyn_mars_violence1D$par$Z)
marsOut$Variable <- c("\\% Supporting two-state solution",
                      "\\% Supporting a one shared state solution",
                      "\\% Supporting a one Islamic state solution",
                      "\\% Saying there is no solution",
                      "\\% Supporting a peace process",
                      "\\% Opposing a peace process",
                      "\\% Supporting military action against Israel",
                      "\\% Opposing military action against Israel",
                      "\\% Supporting suicide bombings",
                      "\\% Opposing suicide bombings",
                      "\\% Very optimistic about peace",
                      "\\% Optimistic about peace",
                      "\\% Pessimistic about peace",
                      "\\% Very pessimistic about peace",
                      "\\% Strongly support Oslo",
                      "\\% Support Oslo",
                      "\\% Oppose Oslo",
                      "\\% Support the \\emph{Intifada}",
                      "\\% Oppose the \\emph{Intifada}",
                      "\\% Who think peace is dead",
                      "\\% Who think the peace process is stalled",
                      "\\% Who think peace is alive",
                      "\\% Support armed attacks generally",
                      "\\% Support armed attacks against civilians",
                      "\\% Support armed attacks against soliders",
                      "\\% Support armed attacks against settlers")


marsOut <- marsOut[  ,c("Variable","V1")]
cat(kable(marsOut, row.names=FALSE, col.names=c("Variable", "Est."),
          caption="ML estimates for latent support for violence",
          digits=2, format="pipe"),
    file="../../Output/Tables/tableD1.txt",
    sep="\n")


## move it into a data frame for merging
gg.extra <- data.frame(ts = 1:dim(dat)[1],
                       month = dat$month,
                       year = dat$year,
                       date = as.yearmon(1994 + seq(0, dim(dat)[1]-1)/12),
                       employStates = c(dyn_mars_uemp$states),
                       emplySE = c(dyn_mars_uemp$states.se),
                       violence1DStates =  c(dyn_mars_violence1D$states[1,]),
                       violence1DSE = c(dyn_mars_violence1D$states.se[1,]))


save(gg.extra,  file="../../Data/ExtraFactors.rdata")

