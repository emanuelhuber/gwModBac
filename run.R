#!/usr/bin/env Rscript
# rm(list=ls(all=TRUE))

##--- ARGUMENTS
args <- commandArgs(trailingOnly <- TRUE)

#args <- c("sim_name", 20, 50, 10, 4, "NULL", 100, 200, 20, "OK")



# arg1 -> outputfile, e.g. "output.txt"
foutput <- "rea"
# arg2 -> nx (= number of cells in x-direction)
# arg3 -> ny (= number of cells in y-direction)
# arg4 -> nz (= number of cells in z-direction)
# arg5 -> number of simulation
# arg6 -> seed
# arg7 -> nx ref
# arg8 -> ny ref
# arg9 -> nz ref
# arg10 -> plot
numOfSim <-  1                    # number of simulations
mySeed <- NULL

DIR <- file.path("/media/huber/Elements/UNIBAS/software/", 
                 "codeR/multilevelMCMC/gwModBac")


if(length(args) >= 1){
  fout <- as.character(args[1])
  initialOptions <- commandArgs(trailingOnly = FALSE)
  fname <- "--file="
  scriptName <- sub(fname, "", initialOptions[grep(fname, initialOptions)])
  DIR <- dirname(scriptName)
}
foutput <- paste0(fout, ".txt")
if(length(args) >= 2){
  nx <- as.integer(args[2])
}
if(length(args) >= 3){
  ny <- as.integer(args[3])
}
if(length(args) >= 4){
  nz <- as.integer(args[4])
}
if(length(args) >= 5){
  numOfSim <- abs(as.integer(args[5]))
}
if(length(args) >= 6){
  if(as.character(args[6]) == "NULL"){
    mySeed <- NULL
  }else{
    mySeed <- as.numeric(args[6])
  }
}
if(length(args) >= 7){
  nxref <- as.integer(args[7])
  if(nxref < 10)  stop("nx ref (arg #7) must be larger than 10!\n")
  if(nxref < nx) stop("nx must be <=  nx ref\n")
}
if(length(args) >= 8){
  nyref <- as.integer(args[8])
  if(nyref < 10)  stop("ny ref (arg #8) must be larger than 10!\n")
  if(nyref < ny) stop("ny must be <=  ny ref\n")
}
if(length(args) >= 9){
  nzref <- as.integer(args[9])
  if(nzref < 2)  stop("nz ref (arg #9) must be larger than 1!\n")
  if(nzref < nz) stop("nx must be <=  nx ref\n")
}
if(length(args) >= 10){
  onePlot <- TRUE
}else{
  onePlot <- FALSE
}


#--- check if DIR exisits
if(!dir.exists(DIR)){
  cat("", file = foutput, append = FALSE)
  stop("This directory (arg #6) does not exist!\n")
}


setwd(DIR)
getwd()


##------ load/install packages
pckgs <- c("zoo", "compiler", "signal", "RandomFields", "akima", "rgdal",
           "raster", "rgeos", "devtools", "Cairo")
inst_pckgs <- pckgs[ !(pckgs %in% installed.packages()[,"Package"]) ]
if(length(inst_pckgs) > 0) install.packages(inst_pckgs)

githubPckgs <- c("GauProMod")
inst_githubPckgs <- githubPckgs[ !(githubPckgs %in% 
                                        installed.packages()[,"Package"]) ]
if(length(inst_githubPckgs) > 0){
  devtools::install_github(paste0("emanuelhuber/", githubPckgs))
}

pkgs_loaded <- lapply(c(pckgs, githubPckgs), require, character.only = TRUE)



invisible(enableJIT(3))

source("fx/RMODFLOW.R")
source("fx/utilityFunctions.R")

pfx <- format(Sys.time(), "%Y_%m_%d_%H%M%S")


#---- directory project ----#
dirProj <- file.path(getwd(), "simulations")
if(!dir.exists(dirProj)) dir.create(path = dirProj)


##----------------------------- MODEL PARAMETERS -----------------------------##
source("parameters.R")   # read all the model parameters

if(exists("nx")){
  modGrid$nx <- nx
}
if(exists("ny")){
  modGrid$ny <- ny
}
if(exists("nz")){
  modGrid$nz <- nz
}

#--- CHECKS
if(modGrid$nx < 10){
  cat("", file = foutput, append = FALSE)
  stop("nx (arg #2) must be larger than 10!\n")
}
if(modGrid$ny < 10){
  cat("", file = foutput, append = FALSE)
  stop("ny (arg #3) must be larger than 10!\n")
}
if(modGrid$nz < 2){
  cat("", file = foutput, append = FALSE)
  stop("nz (arg #4) must be larger than 1!\n")
}
  
  

##-------------------------------- READ DATA ------------------------------##
timeID <- as.character(read.table("data/timeID.txt", sep = ",", header = FALSE,
                                  stringsAsFactors = FALSE))
# OBSERVATIONS
obs <- list()
# river
obs$riv$h <- read.zoo(file = "data/riverStage.txt", sep =",", header = FALSE, 
                      stringsAsFactors =FALSE, format=c("%d.%m.%y"), 
                      FUN = as.POSIXct)
obs$riv$pos <- as.vector(read.table("data/riverStageStation.txt", sep = ",", 
                          header = TRUE, stringsAsFactors = FALSE))
# groundwater
obs$gw$h <- as.matrix(unname(read.table("data/gwHeads.txt", sep = ",", 
                             header = FALSE, stringsAsFactors = FALSE)))
obs$gw$pos <- as.matrix(read.table("data/gwHeadStations.txt", sep = ",", 
                          header = TRUE, stringsAsFactors = FALSE))
weatherFor <- read.zoo(file = "data/weatherFor.txt", sep =",", header = FALSE, 
                      stringsAsFactors =FALSE, format=c("%d.%m.%y"), 
                      FUN = as.POSIXct)
precRef <- read.zoo(file = "data/precRef.txt", sep =",", header = FALSE, 
                      stringsAsFactors =FALSE, format=c("%d.%m.%y"), 
                      FUN = as.POSIXct)


# forecast time steps
timeForc <- (nstp - nstp_pred + 1):nstp
timePast <- 1:(nstp - nstp_pred)
timeIDFor <- timeID[timeForc]
timeIDPast <- timeID[timePast]


##--------------------------- SPATIAL DISCRETISATION -------------------------##
##--- model grid
b <- obs$riv$pos[,"z"] - (-river$slope) * obs$riv$pos[,"y"]
gwMod <- modGrid3D(modGrid, prec = 2, fun = valleyFloor, a = -river$slope, 
                  b = b)
pz_layer1 <-  0.4 # elevation layer 1 = gwMod[[1]] + pz_layer1
gwMod[[1]] <- gwMod[[1]] + pz_layer1
#plot(gwMod[[3]])



##--- river raster
rivPoly <- SpatialPolygons(list(Polygons(list(Polygon(river$perimeter)), 
                           "p1")), 1L)
# rRiv <- rasterize(rivPoly, gwMod[[1]])
rRiv <- rasterizePolygon(rivPoly, gwMod[[1]])
names(rRiv) <- "river"
gwMod <- stackRaster(gwMod, rRiv)
#plot(rRiv)
rm(rRiv)




##--- CHD raster
rCHD <- gwMod[[1]]
rCHD[] <- NA
rCHD[c(1,nrow(rCHD)),] <- 1
rCHD[!is.na(gwMod[["river"]])] <- NA
names(rCHD) <- "CHD"
gwMod <- stackRaster(gwMod, rCHD)
rm(rCHD)
# plot(gwMod[["CHD"]])
# rect(grid$W["min"], grid$L["min"], grid$W["max"], grid$L["max"])

##------------------------------ BOUNDARY CONDITIONS -------------------------##
##--- drinking water extraction Well
wellFrame <- setWells(gwMod, val = wellExt, timeID)
wellFrame[,timeID[1:7]] <- 0  # well off the first 7 days

##--- River
cellsRiv <- cellsFromExtent(gwMod[["river"]], trim(gwMod[["river"]])) 
#--- river bed elevation
xyRiv <- xyFromCell(gwMod[["river"]], cellsRiv)
rivBedAppz <- b + (-river$stepS) * xyRiv[,"y"] +  0.4 - 0.1
rivBedz0 <- rivBedAppz + floor( xyRiv[,"y"] /river$stepL ) * river$stepH 
rivBedz <- rivBedz0 - river$depth

# plot(xyRiv[,"y"], gwMod[[1]][cellsRiv], type="l", asp = 5)
# lines(xyRiv[,"y"],rivBedz0, col="grey", lty=3)
# lines(xyRiv[,"y"], rivBedz, col = "red")
# lines(xyRiv[,"y"], rivBedz + river$minStage, col = "blue")

#--- riverbed conductance
Cr0 <- (10^runif(length(cellsRiv),-3,-1) )
Cr <- Cr0 * res(gwMod)[1]*res(gwMod)[2]/river$bedT
#gwMod[["river"]][cellsRiv] <- rivBedz
riverFrame <- rivGwMod(gwMod, hrel = obs$riv$h, 
                       rivH0z = rivBedz + river$minStage, 
                       rivBedz = rivBedz, 
                       Cr = Cr*cst_mps2mpd, timeID = timeID)

##--- CHD
rowColCHD <- rowColFromCell(gwMod[["CHD"]], 
                            which(!is.na(values(gwMod[["CHD"]]))))
xyCHD <- xyFromCell(gwMod[["CHD"]], which(!is.na(values(gwMod[["CHD"]]))))
idCells <- cellFromRowCol(gwMod, rowColCHD[,"row"], rowColCHD[,"col"])
zl <- as.vector(extract(gwMod[["lay1.top"]], idCells))
ZCHD <- matrix(zl, nrow = length(idCells), ncol = length(timeID), byrow = FALSE)
zlbot <- as.vector(extract(gwMod[[paste0("lay",nlay(gwMod),".bot")]], xyCHD))
ZCHDbot <- matrix(zlbot, nrow=nrow(xyCHD),ncol=nstp,byrow=FALSE)


##-------------------------------- OBSERVATIONS ------------------------------##
##--- River observations
rowColRiv <- rowColFromCell(gwMod[["river"]], cellsRiv)
colnr <- max(rowColRiv[,2])
subSplRiv <- round(seq(1, to = nrow(gwMod[["river"]]),length.out = nrivObs))
rivVal <- t(riverFrame[riverFrame$col == colnr, timeID][subSplRiv,])
xyRivVal <-  xyFromCell(gwMod[["river"]], cellFromRowCol(gwMod[["river"]], 
                        rownr = subSplRiv, colnr = colnr))
xyRivVal[,"x"] <- xyRivVal[,"x"] + res(gwMod)[1]/2

# Boundary condition derivative for GP simulation
bc <- list(x = cbind(rep(xmax(gwMod),nbc) - res(gwMod)[1],
                    seq(from = ymin(gwMod) + res(gwMod)[2], 
                        to   = ymax(gwMod) - res(gwMod)[2],
                        length.out = nbc)),
           v = cbind(rep(1,nbc),
                    rep(0,nbc)),
           y =  rep(0,nbc),
           sigma = 0)

##-------------- PARTICLES - Bacteria infiltration into aquifer --------------##
xyExt <- xyFromCell(gwMod[[1]], cellFromXY(gwMod[[1]], 
                    xy = c(wellExt$x, wellExt$y)))
wellExt$x <- xyExt[1,1]
wellExt$y <- xyExt[1,2]
xyzExtWellFor <- zylParticles(co = xyExt, ro = 0.20, 
                              n=wellExtPart$nxy, zbot=wellExt$zbot, 
                              ztop=wellExt$ztop, npz=wellExtPart$nz)
partFor <- vector(mod="list",length=length(wellExtPart$t))
for(i in seq_along(wellExtPart$t)){
  partFor[[i]] <- setParticles(gwMod, xyzExtWellFor, 
                                releaseTime=wellExtPart$t[i])
}

##----------------------- MONTE CARLO SAMPLING -------------------------------##
Cw <- matrix(0, nrow = numOfSim, ncol = length(timeForc))
it <- 0                            # start at iteration it

while(it < numOfSim){
  it <- it + 1
  cat("***** SIMULATION N0", it, "******\n")
  cat(format(Sys.time(), "   %Y/%m/%d %H:%M:%S \n"))
  
  idRea <- paste0(fout, "_", sprintf("%04d", it))
  dirRun <- file.path(dirProj, idRea)
  suppressWarnings(dir.create(path = dirRun))
  if(!is.null(mySeed)){
    set.seed(mySeed)
  }
  
  ##------------------- HYDRAULIC PROPERTIES SIMULATION -----------------------#
  ##---- hyd. properties
  # riverbed conductance m2/s
  Cr0    <- 10^rnorm(1, mean = prior$Cr$mean, sd = prior$Cr$sd)
  Cr     <- Cr0 * res(gwMod)[1]*res(gwMod)[2]/river$bedT
  ss     <- runif(1, prior$ss$min, prior$ss$max)        # specific storage
  sy     <- runif(1, prior$sy$min, prior$sy$max)        # specific yield
  poros  <- runif(1, prior$p$min, prior$p$max)          # porosity
  Kvani  <- runif(1, prior$K_vani$min, prior$K_vani$max) 
  K_mean <- 10^rnorm(1, mean = prior$K_mean$mean, sd = prior$K_mean$sd)  
  K_sd   <- runif(1, prior$K_sd$min, prior$K_sd$max)     # m/s
  ##--- covariance model
  # horizontal anisotropy angle
  K_hani <- runif(1, prior$K_hani$min, prior$K_hani$max)         
  # streching ratio horizontal
  K_hstr <- runif(1, prior$K_hstr$min, prior$K_hstr$max)       
  # streching ratio vertical
  K_vstr <- runif(1, prior$K_vstr$min, prior$K_vstr$max)         
  # smoothness parameter Matern cov
  K_nu   <- runif(1, prior$K_nu$min, prior$K_nu$max)               
  # correlation length
  K_l    <- runif(1, prior$K_l$min, prior$K_l$max)            
  # nugget sd
  K_nug  <- runif(1, prior$K_nug$min, prior$K_nug$max)
  ##--- simulation
  #--- hydraulic conductivity
  if(!is.null(mySeed)){
    RFoptions(seed = mySeed)
    modGridRef <- modGrid
    modGridRef$nx <- nxref
    modGridRef$ny <- nyref
    modGridRef$nz <- nzref
    gwModRef <- modGrid3D(modGridRef, prec = 2, fun = valleyFloor, 
                          a = -river$slope, b = b)
    gwModRef[[1]] <- gwModRef[[1]] + pz_layer1
    gwMod <- suppressMessages(suppressWarnings(
                  gpHKgwModMultiScale(gwMod, K_hani, K_hstr, 
                                      K_vstr, K_nu, K_sd,
                                      K_l, K_nug, K_mean, 
                                      cst_mps2mpd, gwModRef )))
  }else{
    gwMod <- suppressMessages(suppressWarnings(gpHKgwMod(gwMod, K_hani,
                                                        K_hstr,
                                                       K_vstr, K_nu, K_sd,
                                                       K_l, K_nug, K_mean, 
                                                       cst_mps2mpd )))
  }
  #--- porosity
  gwMod <- porosity(gwMod, poros)    # porosity
  #--- Zonation for "ss" and "sy"
  att.table <- as.data.frame(list(ID = 1, 
                                ss = ss,   # specific storage
                                sy = sy))  # specific yield
  gwMod <- zonation(gwMod, att.table)
  # river bed conductivity
  riverFrame[,"cond"] <- Cr * cst_mps2mpd
  ##--------------------------------------------------------------------------##

  ##------------------- BOUNDARY CONDITION SIMULATION ------------------------##
  ##--- Forcasting river stage and groundwater heads (10-day ahead)
  hfor <- forecastModel(rivh = as.numeric(riverFrame[1, timeIDPast]),
                        prec = as.numeric(precRef),
                        gwh  = obs$gw$h[, timePast],
                        weatherFor = weatherFor, convMod = convMod)
 
  riverFrame[, timeIDFor] <- matrix(hfor$riv[timeForc], 
                                    nrow = nrow(riverFrame),
                                    ncol = length(timeForc), byrow = TRUE) -
                            (riverFrame[1, timeID[1]] - riverFrame[, timeID[1]])
  ##--- Specified head boundary conditions simulation (Gaussian Process)
  rivVal <- t(riverFrame[riverFrame$col == colnr, timeID][subSplRiv,])
  hobs   <- list("x" = rbind(obs$gw$pos[,1:2], xyRivVal), 
                 "y" = c(as.vector(t(hfor$gw)), as.vector(unlist(rivVal))),
                 "t" = seq_along(timeID))
  #---- boundary: CHD ----#
  tryAgain <- TRUE
  while(tryAgain){
    h_sig <- runif(1, prior$h_sig$min, prior$h_sig$max)    # (m) -> sd
    covModel <- list(kernel = "matern",        
                     l      = runif(1, prior$h_lx$min, prior$h_lx$max),
                     v      = runif(1, prior$h_vx$min, prior$h_vx$max),
                     h      = runif(1, prior$h_hx$min, prior$h_hx$max),
                     scale  = c(2,1))
    covModelTime <- list(kernel = "gaussian",
                         l      = runif(1, prior$h_lt$min, prior$h_lt$max),
                         h      = runif(1, prior$h_ht$min, prior$h_ht$max))
    covModels <- list(pos = covModel, time = covModelTime)
    GPCHD <- suppressWarnings(gpCond(obs = hobs, targ = list("x" = xyCHD), 
                              covModels = covModels,
                              sigma = h_sig, op = 2, bc = bc ))
    L <- cholfac(GPCHD$cov)
    hCHD <- suppressWarnings(gpSim(GPCHD , L = L))
    colnames(hCHD) <- c("x","y","t","value")
    valCHD <- matrix(hCHD[,"value"], nrow=nrow(rowColCHD), ncol=nstp, 
                    byrow=TRUE)
    testCHD <- sum(valCHD >= ZCHD) + 
                  sum(valCHD <= ZCHDbot)
    cat("   ")
    if(testCHD > 0){
      cat("+ ")
      next
    }else{
      #---STARTING HEADS
      hinit <- akima::interp(x = xyCHD[,1], y= xyCHD[,2], z = valCHD[,1], 
                        xo = xaxis(gwMod), yo = yaxis(gwMod), linear = TRUE)
      rStrH <- gwMod[[1]]
      rStrH[] <- as.vector(hinit$z)
      rStrH[!is.na(gwMod[["river"]])] <- riverFrame[,timeID[1]]
      if(any((gwMod[["lay1.top"]][] - rStrH[]) < 0) ){
        cat("  +  ")
        next
      }else{
        tryAgain <- FALSE
      }
    }
  }
  CHDFrame <- corCHD(gwMod, rowColCHD, valCHD, timeID)
  gwMod <- initialHeads(gwMod, rStrH)
  ##--------------------------------------------------------------------------##
  
  ##--------------------------- MODFLOW SIMULATION ---------------------------##
  idMF <- "simBC"     # ID MODFLOW model 
  wetting <- c("wetfct" = 0.1 , "iwetit" = 5  , "ihdwet" = 0, "wetdry"= 0.8)
  arguments <- list(rs.model = gwMod, 
                       well = wellFrame, 
                       river = riverFrame,
                         chd = CHDFrame,
                          id = idMF, 
                     dir.run = dirRun, 
                   ss.perlen = 5L, 
           tr.stress.periods = as.Date(timeID[-1], timeFormat),
                     wetting = wetting,
              is.convertible = TRUE,
                         uni = uni,
                  timeFormat = timeFormat)
  suppressMessages(do.call(WriteModflowInputFiles, arguments))
  cat("Run MODFLOW...")
  A <- runModflowUsg(dirpath = dirRun, id = idMF, exe = "mfusg")
  
  #--- check!
  if(!any(grepl("normal termination", A, ignore.case=TRUE))){
    cat("MODFLOW failed!!\n")
    it <- it-1
    unlink(dirRun, recursive=TRUE, force=TRUE)
    if(!is.null(mySeed)){
      mySeed <- mySeed + 1
    }
    next
  }
  cat("... OK!\n")
  ##..........................................................................##

  #   heads.info(fHeads)

  ##.............. MODPATH SIMULATION - BACTERIA INFILTRATION ................##
  idMP2 <- "bacteria"
  suppressMessages(writeModpathInputFiles( id             = idMP2,
                         dir.run         = dirRun,
                         optionFlags     = c("StopOption"          = 2, 
                                             "WeakSinkOption"      = 1,
                                             "WeakSourceOption"    = 2,
                                             "TrackingDirection"   = 2,
                                             "ReferenceTimeOption" = 1), 
                         budgetFaceLabel = NULL, 
                         fbud            = file.path(dirRun, 
                                                     paste0(idMF,".bud")),
                         rs.model        = gwMod,
                         particles       = partFor,
                         ReferenceTime   = nstp,
                         unconfined      = TRUE,
                         verbose         = FALSE))
  cat("   Run MODPATH...")
  B <- runModpath(dirpath = dirRun, id = idMP2, exe = "mp6", 
                  batFile = "runModpath.bat")
  
  if(is.null(B) || !any(grepl("normal termination", B, ignore.case=TRUE))){
    cat("MODPATH (2) failed!!\n")
    it <- it-1
    unlink(dirRun, recursive=TRUE, force=TRUE)
    if(!is.null(mySeed)){
      mySeed <- mySeed + 1
    }
    next
  }
#   ext <- extent3D(gwMod)
#   Pend <- readParticles(idMP2, dirRun, r = gwMOD, type="end")
#   Ppath <- readParticles(idMP2, dirRun, r = gwMOD, type="path")
  
  Pend <- readParticles(idMP2, dirRun, type="end")
  Ppath <- readParticles(idMP2, dirRun,  type="path")
 
  if(all(Pend[,"iLay"] == Pend[,"fLay"] &&
         Pend[,"iRow"] == Pend[,"fRow"] &&
         Pend[,"iCol"] == Pend[,"fCol"])){
    cat("MODPATH (2) > particles did not move!!\n")
    it <- it-1
    unlink(dirRun, recursive=TRUE, force=TRUE)
    if(!is.null(mySeed)){
      mySeed <- mySeed + 1
    }
    next
  }
  cat("... OK!\n")
  ##..........................................................................##

  ##.............. SIMULATION MICROBES CONCENTRATION IN WELL .................##
  Cw[it, ] <- microbesSim(r = gwMod[["river"]], Pend = Pend,
                          Ppath = Ppath, rivh = hfor$riv, 
                          a = bactConcPara$a, b = bactConcPara$b, d = 100, 
                          span = 0.4, lambda = bactConcPara$lambda)
  ##..........................................................................##
  
  if(isTRUE(onePlot)){
    fHeads <- file.path(dirRun , paste0(idMF , ".hds"))
    rHeads <- get.heads(fHeads, kper = 1:nstp, kstp = 1, r = gwMod[[1]])
    hi <- paste0("lay1.head.1.", length(timeID))
    
    nt <- length(hfor$riv)
    timeC <- unique(nt + 1 - Pend[,"iTime"])
    # 1. identify particles coming from the river
    test <- extract(gwMod[["river"]], Pend[,c("x","y")], method = "bilinear")
    vtest <- (!is.na(test) & !is.na(Pend[,c("z")]))
  
    ids <- unique(Ppath[,"id"])
    if(length(ids) > 500){
      ids <- ids[sample.int(n = length(ids), size = 500)]
    }
    polyRiv <- rasterToPolygons(gwMod[["river"]], dissolve=TRUE)
    png(filename = paste0(dirRun, ".png"), 
            width = 480, height = 860, pointsize = 12)
      plot(rHeads[[hi]])
      contour(rHeads[[hi]], levels=seq(0,30,by=0.05), add=TRUE)
      sp::plot(polyRiv, col = rgb(145/255, 238/255, 1, 150/255), add = TRUE)
      plotPathXY(Ppath, id = ids)
      points(Pend[vtest,c("x","y")] ,pch=20,col="black", cex = 1.5) 
      points(Pend[vtest,c("x","y")] ,pch=20,col="green", cex = 0.5) 
      points(Pend[!vtest,c("x","y")] ,pch=20,col="red", cex = 1)
    points(wellExt["x"],wellExt["y"], pch=22, bg="blue",cex=1)
      title(tail(timeID, 1))
    dev.off()
  
  }
  cat("   Concentration = ", 
      format(Cw[it, length(timeForc)], 12, scientific = TRUE), 
      "\n")
  unlink(dirRun, recursive=TRUE, force=TRUE)
}



write.table(format(Cw[, length(timeForc)], 12, scientific = TRUE), 
            file = foutput, append = FALSE, quote = FALSE, sep = "\t", 
            col.names = FALSE, row.names = FALSE)

#   plot(gwMod[["river"]])
#   plotPathXY(Ppath[Ppath[,"id"] %in% 
#                       Pend[which(!is.na(test2)),"id"],])
#   plot(gwMod[["river"]])
#   plotPathXY(Ppath[Ppath[,"id"] %in% 
#                       Pend[which(!is.na(test)),"id"],])
#   
#   
#   plotPathXY(Ppath[which(!is.na(test2)),])
# 
# 
#  
#  
#    fHeads <- file.path(dirRun , paste0(idMF , ".hds"))
#   rHeads <- get.heads(fHeads, kper = 1:nstp, kstp = 1, r = gwMod[[1]])
#   CairoPNG(filename = file.path(dirProj, 
#             paste0(idRea,"_", length(timeID), "_.png")), 
#             width = 480, height = 860, pointsize = 12)
#     hi <- paste0("lay2.head.1.", length(timeID))
#     plot(rHeads[[hi]])
#     contour(rHeads[[hi]], levels=seq(0,30,by=0.05), add=TRUE)
#     points(Pend[,c("x","y")],pch=20,col="blue")  # end (final)
#     plotPathXY(Ppath)
#     plotPathXY(Ppath)
#     points(wellExt["x"],wellExt["y"], pch=22, bg="green",cex=1)
#     title(tail(timeID, 1))
#   dev.off()
#   
#   plot(gwMod[["river"]])
