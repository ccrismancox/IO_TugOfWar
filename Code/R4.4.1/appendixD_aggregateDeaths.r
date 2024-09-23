library(data.table)
library(zoo)

pal2000 <-  fread("../../Data/palestinian_deaths_2000_2008.csv")
pal2009 <-  fread("../../Data/palestinian_deaths_2008_2020.csv")

palDeaths <- merge(pal2000, pal2009, by=colnames(pal2000), all=TRUE)
setnames(palDeaths, "Date of death", "date")
palDeaths <- palDeaths[, "date"]

palDeaths[, date := as.yearmon(as.Date(date))]
palDeaths <- palDeaths[, .N, by=date]
setnames(palDeaths, "N", "deaths")
save(palDeaths, file="PalestinianDeaths.rdata")
