#!/usr/bin/env Rscript

# rm(list=ls(all=TRUE))

args <- commandArgs(trailingOnly <- TRUE)

### ADAPT THIS FILE PATH!!!
ROOT   <- "/home/huber/WORK/UNIBAS/"
DIR0   <- "RESEARCH/GW_uncertainty/auswertung/multiLevelMCMC/"
DIR   <- file.path(ROOT,DIR0)
DIR <- "/media/data/huber/Documents/WORK_NEWNEW/multiLevelMCMC"
setwd(DIR)
getwd()


require(chron)
require(Cairo)
require(zoo)
require(plot3D)
require(MASS)
require(RConics)
require(compiler)
require(signal)
library(RConics)
require(RandomFields)
enableJIT(3)

library(devtools)
devtools::install_github("emanuelhuber/GauProMod")
library(GauProMod)

source("RMODFLOW.R")
source("utilityFunctions.R")

pfx <- format(Sys.time(), "%Y_%m_%d_%H%M%S")



##----------------------------- MODEL PARAMETERS -----------------------------##
source("hyperPara_runModGP_02.R")   # read all the model parameters

##-------------------------------- READ DATA ------------------------------##
# forecast time steps
timeForc <- (nstp - nstp_pred + 1):nstp
timePast <- 1:(nstp - nstp_pred)
dirRef <- file.path(getwd(), "refModGP_a", "001")
#--- timeID
timeIDRef <- readRDS(file.path(dirRef, "timeID.rds"))
timeID <- timeIDRef[1:nstp]
timeIDFor <- timeID[timeForc]
timeIDPast <- timeID[timePast]
#--- gw mod
fObsGWHeads <- file.path(dirRef, "obsGWHeads.rds")
if(file.exists(fObsGWHeads)){
  piez$y <- readRDS(fObsGWHeads)
}else{
  gwModRef <- readRDS(file.path(dirRef, "gwMod.rds"))
  fHeadsRef <- file.path(dirRef , paste0("ref3Dss" , ".hds"))
  gwHeadsRef <- get.heads(fHeadsRef, kper = seq_along(timeID), kstp = 1, 
                          r = gwModRef[[1]])
  piez$y <- headts(gwHeadsRef, piez, gwModRef)
  saveRDS(piez$y, file = file.path(dirRef, "obsGWHeads.rds"))
}
#--- obs
obs <- readRDS(file.path(dirRef , "obs.rds"))
# weather forecast
weatherFor <- readRDS(file.path(dirRef, "weatherForecast.rds"))
# precipitations
precRef0 <- readRDS(file.path(dirRef, "precipitations.rds"))
precRef <- precRef0[seq_along(timeID)]

##---------------------------------------------------------------#
## save

save(timeID, file = "timeID.txt", ascii = TRUE)


write.table(t(timeID), file = "data/timeID.txt", sep = ",", row.names = FALSE,
            col.names = FALSE, quote = TRUE)
write.table( piez$y, file = "data/gwHeads.txt", sep = ",", row.names = FALSE,
            col.names = FALSE, quote = FALSE)
write.table((piez$x), file = "data/gwHeadStations.txt", sep = ",", 
            row.names = FALSE, col.names = c("x", "y", "z"), quote = FALSE )
write.table((obs$riv$h), file = "data/riverStage.txt", sep = ",", 
            col.names = FALSE, quote = FALSE )
write.table((obs$riv$pos), file = "data/riverStageStation.txt", sep = ",", 
            row.names = FALSE, quote = FALSE )
weatherForCB <- cbind(weatherFor[[1]], weatherFor[[2]])
write.table((weatherForCB), file = "data/weatherFor.txt", sep = ",", 
            col.names = FALSE, quote = FALSE, row.names = time(weatherForCB))
write.table((precRef), file = "data/precRef.txt", sep = ",", 
            col.names = FALSE, quote = FALSE )           
## READ
timeID <- as.character(read.table("data/timeID.txt", sep = ",", header = FALSE,
                                  stringsAsFactors = FALSE))
# OBSERVATIONS
obs <- list()
# river
obs$riv$h <- read.zoo(file = "data/riverStage.txt", sep =",", header = FALSE, 
                      stringsAsFactors =FALSE, format=c("%d.%m.%y"), 
                      FUN = as.POSIXct)
obs$riv$pos <- read.table("data/riverStageStation.txt", sep = ",", 
                          header = TRUE, stringsAsFactors = FALSE)
# groundwater
obs$gw$h <- as.matrix(unname(read.table("data/gwHeads.txt", sep = ",", 
                             header = FALSE, stringsAsFactors = FALSE)))
obs$gw$pos <- read.table("data/gwHeadStations.txt", sep = ",", 
                          header = TRUE, stringsAsFactors = FALSE)
weatherFor <- read.zoo(file = "data/weatherFor.txt", sep =",", header = FALSE, 
                      stringsAsFactors =FALSE, format=c("%d.%m.%y"), 
                      FUN = as.POSIXct)
precRef <- read.zoo(file = "data/precRef.txt", sep =",", header = FALSE, 
                      stringsAsFactors =FALSE, format=c("%d.%m.%y"), 
                      FUN = as.POSIXct)
                                  
                                  
                                  