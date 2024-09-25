library(sf)
library(zoo)
library(raster)
library(exactextractr)
library(data.table)
library(ncdf4)
rm(list=ls())

palestine1 <- st_read("../../Data/gadm41_PSE.gpkg", layer = "ADM_ADM_1")
palestine0 <- st_read("../../Data/gadm41_PSE.gpkg", layer = "ADM_ADM_0")


gaza <- subset(palestine1, NAME_1=="Gaza")
westbank <- subset(palestine1, NAME_1=="West Bank") 

 
## Using the URL here because the file is 1.2 Gigs. If you have trouble getting it we can share directly.
options(timeout = max(500, getOption("timeout")))
url <- "https://www.psl.noaa.gov/thredds/fileServer/Datasets/gpcc/full_v2020/precip.mon.total.0.25x0.25.v2020.nc"
download.file(url, destfile="./raindata.nc")
rain.brick <- brick("raindata.nc")

datMonth <- data.table()
proj.match <- rep(FALSE, 300)
dates <- seq.Date(as.Date("1993-01-01"), as.Date("2018-12-01"), by="month")
pb <- txtProgressBar(min = 0, max = length(dates), initial = 0) 
for(i in 1:300){
    date <- as.Date(dates[i])
    idx <- paste0("X", gsub("-", ".", as.character(date)))
    rast <- rain.brick[[idx]]
    rast <- rotate(rast)
    proj.match[i] <-  (st_crs(palestine0)== st_crs(rast))
    rainData.westbank <- exact_extract(rast, westbank, 'mean')
    rainData.gaza <- exact_extract(rast, gaza, 'mean')
    rainData.combined <- exact_extract(rast, palestine0, 'mean')
    datOut <- data.table(date= as.yearmon(date),
                         rainfall= rainData.combined,
                         rainfall.gaza = rainData.gaza,
                         rainfall.wb = rainData.westbank)
    datMonth <- rbind(datMonth, datOut)
    setTxtProgressBar(pb,i)
}
close(pb)
cat("Do all the projections match?", all(proj.match), "\n")


datMonth[,month:=months(date)]
datMonth[, dev.gaza :=(mean(rainfall.gaza,na.rm=T)- rainfall.gaza)/(sd(rainfall.gaza)), by=month] 
datMonth[, dev.wb := (mean(rainfall.wb,na.rm=T)- rainfall.wb)/sd(rainfall.wb), by=month] 
datMonth[, dev := (mean(rainfall,na.rm=T)- rainfall)/sd(rainfall), by=month] 

save(datMonth, file="../../Data/rainData.rdata")
system("rm raindata.nc")
