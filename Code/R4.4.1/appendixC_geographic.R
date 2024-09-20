#' ---
#' author: Casey Crisman-Cox and Michael Gibilisco
#' date; Sept 2024
#' title: Survey descriptions
#' ---
library(ggplot2)
rm(list=ls())

## LOAD DATA
jmcc0 <- subset(read.csv("../../Data/jmcc.csv"), year <= 2018)
jmccWB <- subset(read.csv("../../Data/jmcc_WB.csv"), year <= 2018)
jmccGZ <- subset(read.csv("../../Data/jmcc_GAZA.csv"), year <= 2018)
jmccNew <- subset(read.csv("../../Data/jmcc_2023.csv"), year <= 2018)

cpsrWB <- subset(read.csv("../../Data/cpsr_WB.csv"), year <= 2018)
cpsrGZ <- subset(read.csv("../../Data/cpsr_GAZA.csv"), year <= 2018)
cpsr0 <- subset(read.csv("../../Data/cpsr.csv"), year <= 2018)


## SAMPLE SIZES: descriptive statistics for JMCC
tab <- rbind(c(summary(jmccNew$N[!(is.na(jmccNew$trust_1) & is.na(jmccNew$legis_1))])),
             summary(jmccWB$N[!(is.na(jmccNew$trust_1) & is.na(jmccNew$legis_1))]), 
             summary(jmccGZ$N[!(is.na(jmccNew$trust_1) & is.na(jmccNew$legis_1))]))
tab <- round(tab, digits=1)
rownames(tab) <- c("JMCC:All","JMCC:WB", "JMCC:GZ")
print(tab)
## proportion from west bank on average

print(round(mean(jmccWB$N[!(is.na(jmccNew$trust_1) & is.na(jmccNew$legis_1))]/
                          jmccNew$N[!(is.na(jmccNew$trust_1) & is.na(jmccNew$legis_1))])*100))



## SAMPLE SIZES: descriptive statistics for PCPSR
tab <- rbind(c(summary(cpsr0$size[!is.na(cpsr0$supportHamas)])),
             summary(cpsrWB$size[!is.na(cpsr0$supportHamas)]), 
             summary(cpsrGZ$size[!is.na(cpsr0$supportHamas)]))
tab <- round(tab, digits=1)
rownames(tab) <- c("PCPSR:All","PCPSR:WB", "PCPSR:GZ")
tab
summary(cpsrWB$size[!is.na(cpsr0$supportHamas)]/cpsr0$size[!is.na(cpsr0$supportHamas)]*100)



## Graphs

myx <- as.Date(paste0(jmcc0$year,"/",jmcc0$month,"/",1),"%Y/%B/%d")
#myx <- 1:dim(jmcc0)[1]
ggdat <- data.frame(y = c(jmcc0$trust_1, jmccWB$trust_1, jmccGZ$trust_1,
                          jmcc0$trust_2, jmccWB$trust_2, jmccGZ$trust_2,
                          cpsr0$supportFatah, cpsrWB$supportFatah, cpsrGZ$supportFatah,
                          cpsr0$supportHamas, cpsrWB$supportHamas, cpsrGZ$supportHamas,
                          jmcc0$legis_1, jmccWB$legis_1, jmccGZ$legis_1,
                          jmcc0$legis_2, jmccWB$legis_2, jmccGZ$legis_2),
                    x = rep(myx, 18),
                    Area = rep(rep(c("Both", "West Bank", "Gaza Strip"), each=dim(jmcc0)[1]),6),
                    Actor = rep(c("Fatah", "Fatah", "Fatah", 
                                  "Hamas", "Hamas", "Hamas", 
                                  "Fatah", "Fatah", "Fatah", 
                                  "Hamas", "Hamas","Hamas",
                                  "Fatah", "Fatah", "Fatah", 
                                  "Hamas", "Hamas","Hamas"), each=dim(jmcc0)[1]),
                    measure = rep(rep(c("Trust","Support","Vote"),each=6), each = dim(jmcc0)[1]))
ggdat <- ggdat[!is.na(ggdat$y),]
ggdat$measure <- factor(ggdat$measure, levels=c("Trust","Support","Vote"), ordered=T)

ggano <- data.frame(val = c(
                   paste("WB-Both cor: ", round(cor(jmcc0$trust_1, jmccWB$trust_1, use="complete.obs"), digits=2)),
                   paste("GZ-Both cor: ",round(cor(jmcc0$trust_1, jmccGZ$trust_1, use="complete.obs"), digits=2)),
                   paste("WB-Both cor: ",round(cor(jmcc0$trust_2, jmccWB$trust_2, use="complete.obs"), digits=2)),
                   paste("GZ-Both cor: ",round(cor(jmcc0$trust_2, jmccGZ$trust_2, use="complete.obs"), digits=2)),
                   paste("WB-Both cor: ",round(cor(cpsr0$supportFatah, cpsrWB$supportFatah, use="complete.obs"), digits=2)),
                   paste("GZ-Both cor: ",round(cor(cpsr0$supportFatah, cpsrGZ$supportFatah, use="complete.obs"), digits=2)),
                   paste("WB-Both cor: ",round(cor(cpsr0$supportHamas, cpsrWB$supportHamas, use="complete.obs"), digits=2)),
                   paste("GZ-Both cor: ",round(cor(cpsr0$supportHamas, cpsrGZ$supportHamas, use="complete.obs"), digits=2)),
                   paste("WB-Both cor: ",round(cor(jmcc0$legis_1, jmccWB$legis_1, use="complete.obs"), digits=2)),
                   paste("GZ-Both cor: ",round(cor(jmcc0$legis_1, jmccGZ$legis_1, use="complete.obs"), digits=2)),
                   paste("WB-Both cor: ",round(cor(jmcc0$legis_2, jmccWB$legis_2, use="complete.obs"), digits=2)),
                   paste("GZ-Both cor: ",round(cor(jmcc0$legis_2, jmccGZ$legis_2, use="complete.obs"), digits=2))),
          x = rep(as.Date("1998-02-01"), 12),
          y = c(16,11,56,51,16,11,56,51,16,11,16,11),
          Actor=rep(c("Fatah","Fatah","Hamas", "Hamas"), 3),
          measure=rep(c("Trust","Support","Vote"),each=4),
          Area=NA)
ggano$measure <- factor(ggano$measure, levels=c("Trust","Support","Vote"), ordered=T)




fig <- ggplot(ggdat, aes(x=x, y=y, color=Area,shape=Area)) + 
  geom_line(linewidth=1.1) + geom_point(size=2) + facet_grid(measure~Actor) + theme_bw(18) + 
  scale_color_manual(values=c("#fc8d62", "#66c2a5",  "#8da0cb")) +
  scale_x_date(breaks=as.Date(paste(seq(from=1994,to=2018, by=4), "-02-01",sep="")),
               labels = seq(from=1994,to=2018, by=4)) +
  ylab("Percent respondents") + xlab("Date") + 
  geom_text(data=ggano, aes(x=x,y=y,label=val),color="Black",size=3.5)+
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = "bottom",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-20,-20,-9,-20))

ggsave("../../Output/Figures/figureC2.pdf", plot = fig, 
       width = 7.75, height = 7, units = "in")

  

  
  

