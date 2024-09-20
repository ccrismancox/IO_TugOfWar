#' ---
#' author: Casey Crisman-Cox and Michael Gibilisco
#' date: Sept. 2024
#' title: civilian and non-civilian attacks
#' ---

rm(list=ls())

load("../../Data/actionsSetup_byAttackType.Rdata")
regData$lag.Hattacks.mil <- (dat$lag.Hattacks.mil >0)*1
regData$lag.Fattacks.mil <- (dat$lag.Fattacks.mil >0)*1
regData$lag.Hattacks.notciv <- (dat$lag.Hattacks.notciv >0)*1
regData$lag.Fattacks.notciv <- (dat$lag.Fattacks.notciv >0)*1
regData$lag.Hattacks.civ <- (dat$lag.Hattacks.civ >0)*1
regData$lag.Fattacks.civ  <- (dat$lag.Fattacks.civ >0)*1

regData$Hsoft <- (dat$Hsoft >0)*1
regData$Fsoft  <- (dat$Fsoft >0)*1

sum(regData$lag.Fattacks==1)
sum(regData$lag.Fattacks.civ==1)
sum(regData$Fsoft==1)
sum(regData$lag.Fattacks.notciv==1)
(23-13)/23 #percentage change in moving to only civ sample
(23-15)/23



sum(regData$lag.Hattacks==1)
sum(regData$lag.Hattacks.civ==1)
sum(regData$Hsoft==1)
sum(regData$lag.Hattacks.notciv==1)
(131-99)/131
(131-82)/131


sum(dat$Hattacks.civ)/sum(dat$Hattacks)
sum(dat$Fattacks.civ)/sum(dat$Fattacks)

mean(dat[Hattacks.civ==1]$Hkills.civ)
mean(dat[Fattacks.civ==1]$Fkills.civ)



regData <- merge(regData, regData0, by=colnames(regData)[colnames(regData) %in% colnames(regData0)])
summary(regData)

mod0.notciv <- lm(diff.states~lag.Hattacks.notciv + lag.Fattacks.notciv + second.int + 
                    L.diff.emp + L.diff.violence + timesinceelection + L.diff.states + 
                    lag.Hattacks.notciv:lag.states + lag.Fattacks.notciv:lag.states, data=regData,
           x=T,y=T)
summary(mod0.notciv)


mod0.civ <- lm(diff.states~lag.Hattacks.civ + lag.Fattacks.civ + second.int + 
                    L.diff.emp + L.diff.violence + timesinceelection + L.diff.states + 
                    lag.Hattacks.civ:lag.states + lag.Fattacks.civ:lag.states, data=regData,
               x=T,y=T)
summary(mod0.civ)
regData_actions <- copy(regData)


mod0.soft <- lm(diff.states~Hsoft+ Fsoft + second.int + 
                    L.diff.emp + L.diff.violence + timesinceelection + L.diff.states + 
                    lag.Hattacks.notciv:lag.states + lag.Fattacks.notciv:lag.states, data=regData,
                  x=T,y=T)
summary(mod0.soft)



save(list=c("mod0.civ", "mod0.notciv"), file="../Results_220/subsetactions_robust.rdata")

