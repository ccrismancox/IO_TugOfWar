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

## using the URL here because the file is 1.2 Gigs. If you have trouble getting it we can share directly.
options(timeout = max(500, getOption("timeout")))
url <- "https://www.psl.noaa.gov/thredds/fileServer/Datasets/gpcc/full_v2020/precip.mon.total.0.25x0.25.v2020.nc"
download.file(url, destfile="./raindata.nc")
rain.brick <- brick("raindata.nc")

datMonth <- data.table()

for(i in seq.Date(as.Date("1993-01-01"), as.Date("2018-12-01"), by="month")){
    i <- as.Date(i)
    idx <- paste0("X", gsub("-", ".", as.character(i)))
    rast <- rain.brick[[idx]]
    extent(rast)
    rast <- rotate(rast)
    cat("Matching projections in month", as.character(i), "\t", st_crs(palestine0)== st_crs(rast), "\n")
    rainData.westbank <- exact_extract(rast, westbank, 'mean')
    rainData.gaza <- exact_extract(rast, gaza, 'mean')
    rainData.combined <- exact_extract(rast, palestine0, 'mean')
    datOut <- data.table(date= as.yearmon(i),
                         rainfall= rainData.combined,
                         rainfall.gaza = rainData.gaza,
                         rainfall.wb = rainData.westbank)
    datMonth <- rbind(datMonth, datOut)
}
datMonth[,month:=months(date)]
datMonth[, dev.gaza :=(mean(rainfall.gaza,na.rm=T)- rainfall.gaza)/(sd(rainfall.gaza)), by=month] 
datMonth[, dev.wb := (mean(rainfall.wb,na.rm=T)- rainfall.wb)/sd(rainfall.wb), by=month] 
datMonth[, dev := (mean(rainfall,na.rm=T)- rainfall)/sd(rainfall), by=month] 

save(datMonth, file="../../Data/rainData.rdata")
system("rm raindata.nc")
