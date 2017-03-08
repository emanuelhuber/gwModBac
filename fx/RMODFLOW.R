
# require(akima)
# require(rgdal)
# require(raster)
# require(rgeos)


 #cat(getwd())
#-------------------------------------------------------------# 
#source(file.path("..", "utils/binaryConnections.R"))
##--- file position (same function name as in MATLAB)
#  returns the current position in the specified file/connection con
ftell <- function(con){
  return(seek(con))
}
# Move to specified position in file
# origin : - current, start or end
fseek <- function(con, where, origin="current"){
  seek(con, where=where,origin=origin)
}
# number of bytes in connection
# file.info(filename)$size
flen <- function(con){
  pos0 <- ftell(con)
  seek(con,0,"end")
  pos <- ftell(con)
  seek(con,where=pos0,"start")
  return(pos)
}
# place pointor at the begining
frewind <- function(con){
  seek(con,where="0",origin="start")  # start position
}
# test end of file
feof <- function(con){
  return(ftell(con) == flen(con))
}
#----------------------------------------------------------------#



#----------------------------------------------------------------#
#source(file.path("..", "GIS/polygonizer.R"))
# need to install "gdal_polygonize.py" > install python-gdal
polygonizer <- function(x, outshape=NULL, gdalformat = 'ESRI Shapefile', 
                        pypath=NULL, readpoly=TRUE, quietish=TRUE) {
  # x: an R Raster layer, or the file path to a raster file recognised by GDAL
  # outshape: the path to the output shapefile (if NULL, a temporary file will 
  # be created)
  # gdalformat: the desired OGR vector format
  # pypath: the path to gdal_polygonize.py (if NULL, an attempt will be made to 
  # determine the location
  # readpoly: should the polygon shapefile be read back into R, and returned by 
  # this function? (logical)
  # quietish: should (some) messages be suppressed? (logical)
  if (isTRUE(readpoly)) require(rgdal)
  if (is.null(pypath)) {
    pypath <- Sys.which('gdal_polygonize.py')
  }
  if (!file.exists(pypath)){
    stop("Can't find gdal_polygonize.py on your system.")
  }
  owd <- getwd()
  on.exit(setwd(owd))
  setwd(dirname(pypath))
  if (!is.null(outshape)) {
    outshape <- sub('\\.shp$', '', outshape)
    f.exists <- file.exists(paste(outshape, c('shp', 'shx', 'dbf'), sep='.'))
    if (any(f.exists)) 
      stop(sprintf('File already exists: %s', 
                   toString(paste(outshape, c('shp', 'shx', 'dbf'), 
                                  sep='.')[f.exists])), call.=FALSE)
  } else outshape <- tempfile()
  if (is(x, 'Raster')) {
    require(raster)
    writeRaster(x, {f <- tempfile(fileext='.asc')})
    rastpath <- normalizePath(f)
  } else if (is.character(x)) {
    rastpath <- normalizePath(x)
  } else stop('x must be a file path (character string), or a Raster object.')
  system2('python', args=(sprintf('"%1$s" "%2$s" %3$s-f "%4$s" "%5$s.shp"', 
                                  pypath, rastpath, ifelse(quietish, '-q ', ''),
 
                                  gdalformat, outshape)))
  if (isTRUE(readpoly)) {
    shp <- readOGR(dirname(outshape), layer = basename(outshape), 
                    verbose = !quietish)
    return(shp) 
  }
  return(NULL)
}
#----------------------------------------------------------------#



# - bbox3D()
# - plotPath3D()
# - raster to matrix: format(raster::as.matrix(r), justify="right", width=2)
# - cellSaturation()
#       # A hydraulic head value that exceeds land-surface elevation indicates 
# cell saturation (ﬁg. 18).
#       r <- rs.run [["lay1.head"]] > rs.model [["lay1.top"]] & 
#               rs.heads [[1]] < 1e30
#       r <- ratify (r)
#       levels(r) <- cbind(levels(r)[[1]], 
#                       att = c("partially saturated", "saturated"))
#       names(r) <- "lay1.saturated"
#       rs.run <- stack (rs.run , r)
# - gwTable()
#       # Determine the simulated elevation of the water table (ﬁg. 15). 
#           Heads greater than the land surface elevation
#       # are speciﬁed at land surface.
# x1 = raster stack
# x2 = raster to add/stack to x1
#    if there is already a x2 raster, it will be overwritten

ExportRasterStack <- function (rs, filepath, zip = "", 
                                col = rainbow(250, start = 0, end = 0.8)){
  dir.create(filepath, showWarnings = FALSE, recursive = TRUE)
  dir.create(filepath.csv <- file.path(filepath,  "CSV"), showWarnings = FALSE)
  dir.create(filepath.png <- file.path(filepath,  "PNG"), showWarnings = FALSE)
  dir.create(filepath.tif <- file.path(filepath, "TIFF"), showWarnings = FALSE)
  dir.create(filepath.rda <- file.path(filepath,  "RDA"), showWarnings = FALSE)
  dir.create(filepath.kml <- file.path(filepath,  "KML"), showWarnings = FALSE)
  n <- 0L
  for (i in names(rs)) {
      n <- n + 1L
      fig.num <- formatC(n, width = 2, format = "d", flag = "0")
      f <- file.path(filepath.csv, paste(fig.num, "_", i, ".csv", 
          sep = ""))
      m <- matrix(data = rs[[i]][], nrow = nrow(rs), ncol = ncol(rs), 
          byrow = TRUE)
      write.table(m, file = f, quote = FALSE, sep = ",", na = "", 
          row.names = FALSE, col.names = FALSE, qmethod = "double")
      f <- file.path(filepath.png, paste(fig.num, "_", i, ".png", 
          sep = ""))
      png(filename = f, width = 7, height = 7, units = "in", 
          pointsize = 12, res = 1200, antialias = "subpixel")
      plot(rs[[i]], maxpixels = length(rs[[i]]), col = col, 
          main = names(rs[[i]]), asp = 1)
      dev.off()
      f <- file.path(filepath.tif, paste(fig.num, "_", i, ".tif", 
          sep = ""))
      writeRaster(rs[[i]], filename = f, format = "GTiff", 
          overwrite = TRUE, NAflag = -999)
  }
  base.name <- "raster"
  f <- file.path(filepath.rda, "rasters.rda")
  save(rs, file = f)
  f <- file.path(filepath.kml, "rasters.kml")
  crs <- "+proj=longlat +datum=WGS84"
  rs <- projectRaster(rs, crs = crs, method = "ngb", alignOnly = FALSE)
  KML(rs, f, col = col, maxpixels = ncell(rs) * 2, blur = 5, 
      zip = zip, overwrite = TRUE)
  invisible()
}
      
# to remove a raster, 
stackRaster <- function(r1, r2, verbose=FALSE, ...){
  if(names(r2) %in% names(r1)){
    r1[[names(r2)]] <- r2
    if(isTRUE(verbose)){
      cat("I overwrite", names(r2), "\n")
    }
  }else{
    r1 <- stack(r1,r2,...)
  }
  return(r1)
}

dropRaster <- function(r, i){
  if(is.character(i)){
    if(i %in% names(r)){
      r <- dropLayer(r, i)    
    }else{
      stop(paste0("raster '", i, 
          "' does not exist in r!\n"))
    }
  }else{
     stop("error")  
  }
  return(r)
}

# Removing outer rows and columns that are all inactive results in the 
# horizontal model grid. A summary of the model grid attributes is 
# shown in table 1.
# x = raster stack
# y = raster that define the grid size

trimRaster <- function(x, y){
  modelExtent <- trim(y)
  x <- stack(lapply(names(x), .cropRaster, x, modelExtent), quick = TRUE)
  return(x)
}

### CHECK EXTRACT !!!!!
valueFromRowCol <- function(r, rc){
  xx <- raster::cellFromRowCol(r, rc[,1], rc[,2])
  return(raster::extract(r, xx))
}
valueFromRowColCombine <- function(r, rc){
  xx <- raster::cellFromRowColCombine(r, rc[,1], rc[,2])
  return(raster::extract(r, xx))
}
valueFromXY <- function(r, xy){
  xx <- raster::cellFromXY(r, xy)
  return(raster::extract(r, xx))
}
valueFromY <- function(r, y){
  rowx <- raster::rowFromY(r, y)
  xx <- raster::cellFromRow(r, rowx)
  return(raster::extract(r, xx))
}
valueFromX <- function(r, x){
  colx <- raster::colFromX(r, x)
  xx <- raster::cellFromCol(r, colx)
  return(raster::extract(r, xx))
}
valueFromLine <- function(r, lns){
  xx <- raster::cellFromLine(r, lns)
  return(raster::extract(r, xx))
}
valueFromPolygon <- function(r, p, weights=FALSE){
  xx <- raster::cellFromPolygon(r, p, weights=weights)
  return(raster::extract(r, xx))
}
 

.cropRaster <- function(i, x , modelExtent){
  return(crop(x[[i]], modelExtent))
}

leftJoin <- function(x, y, by){
  cbind(merge(x, y, by=by))
}

plotHeads<-function(x, cont=TRUE,vec=FALSE){
  val <- x[]
  if(all(is.na(val))){
  cat("all values are 'NA', no plot!\n")
  }else{
  plot(x)
  if(cont){
    contour(x, add=TRUE)
  }
  if(vec){
    plotVec(x)
  }
  }
}

nlay <- function(gwMod){
  rasterNames <- names(gwMod)
    lbot <- grep('^lay[[:digit:]]+.bot$',rasterNames) # match layi.bot
    nlay <- length(lbot)
  return(nlay)
}

# r = rasterStack (rs.model)
bbox3D <- function(r){
  nl <- nlay(r)
  minX <- xmin(r)
  maxX <- xmax(r)
  minY <- ymin(r)
  maxY <- ymax(r)
  minZ <- minValue(r[[paste0("lay",nl,".bot")]])
  maxZ <- maxValue(r[["lay1.top"]])
  return(c(minX, maxX,minY, maxY, minZ, maxZ))
}

xaxis <- function(x){
  xFromCol(x, col=1:ncol(x))
}
    
yaxis <- function(x){
  yFromRow(x, row = 1:nrow(x))
}    
  

extent3D <- function(r){
  aa<- extent(r)
  nl <- nlay(r)
  minZ <- minValue(r[[paste0("lay",nl,".bot")]])
    maxZ <- maxValue(r[["lay1.top"]])
    return(c("xmin" = aa[1], "xmax"=aa[2], "ymin"=aa[3], "ymax"=aa[4], 
             "zmin"=minZ, "zmax"=maxZ))
}

zaxis <- function(x,idCell = 1){
  rnames <- c("lay1.top",paste0("lay", 1:nlay(x),".bot"))
  as.vector(cellStats(x[[rnames]],mean))
}

###############################################################################
##-------------------- UTILITY FUNCTIONS (MODFLOW/MODPATH) ------------------##
###############################################################################


##--- string cleaner
purify <- function(desc) {
  return(tolower(gsub("(^ +)|( +$)", "", desc)))
}



# for windows, indicate the full path of the exe file
# for linux: if the executable is in /usr/bin, 
#             just give the name, else the path
writeBatFile <- function(dirpath=NULL,id=NULL, exe="mfusg", 
                          batFile="runModflow.bat", ext=".nam"){
  if(is.null(dirpath)){
  stop("dirpath does not exist!\n")
  }
  batFile <- file.path (dirpath , batFile)
  cmd <- c(paste("cd", shQuote(dirpath)),
  paste(shQuote(exe), shQuote(paste0(id , ext))))
  cat(cmd , file = batFile , sep = "\n")
  return(batFile)
}


# for windows, indicate the full path of the exe file
# for linux: if the executable is in /usr/bin, 
#           just give the name, else the path
runModflowUsg <- function(dirpath=NULL, id=NULL,exe="mfusg", 
                          batFile="runModflow.bat", ext=".nam"){
  file.bat <- writeBatFile(dirpath=dirpath, id=id, exe=exe, 
                           batFile=batFile, ext=ext)
  Sys.chmod(file.bat , mode = "755")
  output <- system(shQuote(file.bat), intern = TRUE)
  return(output)
}

# for windows, indicate the full path of the exe file
# for linux: if the executable is in /usr/bin, 
#             just give the name, else the path
runModpath <- function(dirpath=NULL, id=NULL, exe="mp6", 
                       batFile="runModpath.bat", ext=".mpsim"){
  file.bat <- writeBatFile(dirpath=dirpath, id=id, exe=exe, 
                           batFile=batFile, ext=ext)
  Sys.chmod(file.bat , mode = "755")
  output <- system(shQuote(file.bat), intern = TRUE)
  return(output)
}


skipBinMarker <- function(con,nskip){
  if(nskip > 0) invisible(readBin(con, "integer", n=1L, size=nskip))
}


readEndParticles <- function(f){
  read.table(f, skip=5, header=FALSE)
}



###############################################################################
##------------------------ READ MODFLOW BINARIES: HEADS ---------------------##
###############################################################################

getValidHeadDesc <- function(){
  valid.desc <- c("head", "drawdown", "subsidence", "compaction",
                    "critical head", "head in hgu", "ndsys compaction",
                    "z displacement", "d critical head", "layer compaction",
                    "dsys compaction", "nd critical head", "system compaction",
                    "preconsol stress", "change in pcstrs", "effective stress",
                    "change in eff-st", "void ratio", "thickness",
                    "center elevation", "geostatic stress", "change in g-strs")
  return(valid.desc)
}

bytesOffsetHeads <- function(con){
    valid.desc <- getValidHeadDesc()
  for(nskip in 0:4){
  seek(con,where="0",origin="start")  # start position
  if(nskip > 0) skp <- readBin(con, "integer", n=nskip, size=1L)
#   skipBinMarker(con, nskip)
  kstp  <- readBin(con, "int", n=1L, size=4L)
  kper  <- readBin(con, "int", n=1L, size=4L)
  pertim <- readBin(con, "numeric", n=1L, size=4L)
  totim  <- readBin(con, "numeric", n=1L, size=4L)
  desc <- readBin(readBin(con, "raw", n=16L, size=1L), "character", n=1L)
  desc <- purify(desc)
  if(!desc %in% valid.desc){
      #  cat(paste("variable description not recognized:", desc))
  }else{
#     cat("desc = ",desc, "\n")
#     cat("bytes offset = ",nskip, "\n")
    seek(con,where="0",origin="start")  # start position
    return(as.integer(nskip))
    break
  }
  }
  stop("no bytes offset found!\n")
}  
  
  
# return information for each layers of heads
# (kper, kstep, pertim, totim, desc...
heads.heading <- function(con, nskip){
  skipBinMarker(con, nskip)
#   kskip <- readBin(con, "int", n=1L, size=4L)    # skp
  kstp  <- readBin(con, "int", n=1L, size=4L)
  kper  <- readBin(con, "int", n=1L, size=4L)
  pertim <- readBin(con, "numeric", n=1L, size=4L)
  totim  <- readBin(con, "numeric", n=1L, size=4L)
  desc <- readBin(readBin(con, "raw", n=16L, size=1L), "character", n=1L)
  desc <- purify(desc)
  nc   <- readBin(con, "integer", n=1L, size=4L)  # ncol
  nr   <- readBin(con, "integer", n=1L, size=4L)  # nrow
  nl   <- readBin(con, "integer", n=1L, size=4L)  # nlayer
  return(c(kstp=kstp, kper=kper, pertim=pertim,totim=totim,
       desc=desc, ncol=nc,     nrow=nr,     nlay=nl))
}
# return the heads for the chosen layers
heads.layer <- function(con, nskip){
  skipBinMarker(con, nskip)
#   kskip <- readBin(con, "int", n=1L, size=4L)    # skp
  kstp  <- readBin(con, "int", n=1L, size=4L)
  kper  <- readBin(con, "int", n=1L, size=4L)
  pertim <- readBin(con, "numeric", n=1L, size=4L)
  totim  <- readBin(con, "numeric", n=1L, size=4L)
  desc <- readBin(readBin(con, "raw", n=16L, size=1L), "character", n=1L)
  desc <- purify(desc)
  nc   <- readBin(con, "integer", n=1L, size=4L)  # ncol
  nr   <- readBin(con, "integer", n=1L, size=4L)  # nrow
  nl   <- readBin(con, "integer", n=1L, size=4L)  # nlayer
  skipBinMarker(con, nskip)
  skipBinMarker(con, nskip)
#   skp <- readBin(con, "int", n=1L, size=4L)    # skp
#   skp <- readBin(con, "int", n=1L, size=4L)
  v      <- readBin(con, "numeric", n=nr * nc, size=4L)
    d <- matrix(v, nrow=nr, ncol=nc, byrow=TRUE)
    return(list(d=d, kstp=kstp, kper=kper, desc=desc,
                nlay=nl, pertim=pertim, totim=totim))
}
# user selection: row, col, layer, kper, kstep
# r is a raster corresponding to the head layers
# return - a list of list of matrix head for each layer, kper, kstp
#     - is r is defined, retunr a raster stack
get.heads <- function(f, kper=NULL,kstp=NULL,nl=NULL, r=NULL, rclip=NULL, 
                      hdry =  -888){
  con <- file(f, open="rb", encoding="bytes")
  nskip   <- bytesOffsetHeads(con)
  Nbyt <- flen(con)      # number of bytes
  dframe <- heads.info(f)
  nbyt <- Nbyt/nrow(dframe)   # number of records in file
  KPER <- as.numeric((dframe$kper))
  KSTP <- as.numeric((dframe$kstp))
  NLAY <-as.numeric((dframe$nlay))
  NCOL <- as.numeric((dframe$ncol))
  NROW <- as.numeric((dframe$nrow))
  if(is.null(kper)){
    kper <- max(KPER)
  }else if( !all(kper %in% KPER)){
    close(con, type="rb")
    stop("kper beyond range! Check with heads.info()!\n")
  }
  if(is.null(kstp)){
    kper <- max(KSTP)
  }else if( !all(kstp %in% KSTP)){
    close(con, type="rb")
    stop("kstp beyond range! Check with heads.info()!\n")
  }
  if(is.null(nl)){
    nl <- 1:max(NLAY)
#     kper <- max(NLAY)
  }else if( !all(nl %in% NLAY)){
    close(con, type="rb")
    stop("nl beyond range! Check with heads.info()!\n")
  }
  # selection
  sel <- which((KPER %in% kper) & (KSTP %in% kstp) & (NLAY %in% nl))
  
  seek(con,where="0",origin="start")  # start position

  lst <- vector(mode="list",length=length(sel))
  for(i in seq_along(sel)){
     # % move to start of last layer in file
    seek(con, where=(sel[i]-1)*nbyt,origin="start")  
    lst[[i]] <- heads.layer(con, nskip)
  }
  close(con, type="rb")
  if(!is.null(r) && NCOL==ncol(r) && NROW==nrow(r)){
    rs <- stack()
    r0 <- r
    for(i in seq_along(sel)){
      r0[] <- lst[[i]]$d
      names(r0) <- paste0("lay",lst[[i]]$nlay,".",lst[[i]]$des)
#       if(length(unique(lst[[i]]$kstp)) > 1){
        names(r0) <- paste0(names(r0), ".",lst[[i]]$kstp)
#       }
#       if(length(unique(lst[[i]]$kper)) > 1){
        names(r0) <- paste0(names(r0), ".",lst[[i]]$kper)
#       }
      if(!is.null(rclip) && dim(rclip)[3]==max(NLAY)){
        namesClip <- names(rclip)
        idRClip <- grep(paste("lay",lst[[i]]$nlay,sep=""), namesClip)
        if(length(idRClip)){
          cat("Clipping",names(r0),"\n")
          r0[is.na(rclip[[idRClip]])] <- NA
        }
      }
      r0[r0 == hdry] <- NA
      rs <- stack(rs,r0)
    }
    return(rs)
    
  }else{
    return(lst)
  }
}
heads.info <- function(f){
  con <- file(f, open="rb", encoding="bytes")
  nskip   <- bytesOffsetHeads(con)
#   kskip <- readBin(con, "int", n=1L, size=4L)    # skp
  skipBinMarker(con, nskip)
  kstp <- readBin(con, "int", n=1L, size=4L)
  kper <- readBin(con, "int", n=1L, size=4L)
  nbytes <- 4L
  pertim <- readBin(con, "numeric", n=1L, size=nbytes)
  totim  <- readBin(con, "numeric", n=1L, size=nbytes)
  desc <- readBin(readBin(con, "raw", n=16L, size=1L), "character", n=1L)
  desc <- purify(desc)
  nc   <- readBin(con, "integer", n=1L, size=4L)
  nr   <- readBin(con, "integer", n=1L, size=4L)
  nl   <- readBin(con, "integer", n=1L, size=4L)
  skipBinMarker(con, nskip)
  skipBinMarker(con, nskip)
#   skp <- readBin(con, "int", n=1L, size=4L)    # skp
#   skp <- readBin(con, "int", n=1L, size=4L)    # skp
    
  n <- ftell(con)    # n = 10*4 + 16 bytes
  readBin(con, "numeric", n=1L, size=4L)
  floatlen <- ftell(con)-n
  fseek(con, where= floatlen*(nr * nc * abs(nl) - 1 ))
  kper <- 4L
  skipBinMarker(con, nskip)
#   readBin(con, "int", n=1L, size=4L)    # skp

  # 648232 + 4 bytes

  nbyt <- ftell(con)  # 56 + 648232 +4
  # seek(con,NA,"end")
  Nbyt <- flen(con)    # number of bytes
  seek(con,where="0",origin="start")  # start position

  nRec <- Nbyt/nbyt  # number of records in file

  # read heading for each layers
  lst <- list()
  dframe <- as.data.frame(matrix(nrow=nRec, ncol=8))
  for(i in 1:nRec){
    # % move to start of last layer in file
    seek(con, where=(i-1)*nbyt,origin="start")
    v <- heads.heading(con,nskip)
    dframe[i,] <- v
  }
  names(dframe) <- names(v)
  close(con, type="rb")
  return(dframe)
}

display.heads.info <- function(dframe){
  budlabels=unique(dframe$desc)
  #iLabInFile=NaN(nRec,1);
#   for i=1:nRec, iLabInFile(i)=strmatchi(B(i).label,budlabels); end
  nRec <- nrow(dframe)
  periodsInFile= as.numeric(dframe$kper)
  NPER = length(unique(periodsInFile))

  tstpsInFile  = as.numeric(dframe$kstp)
  NSTP = length(unique(tstpsInFile))

  NLBL = length(budlabels)     # already unique
  NCOL = as.numeric(dframe$ncol[nrow(dframe)])
  NROW = as.numeric(dframe$nrow[nrow(dframe)])
  NLAY = as.numeric(dframe$nlay[nrow(dframe)])

  pif=unique(periodsInFile)
  tif=unique(tstpsInFile)

  cat('File contains the following:\n');
  cat('Stress Periods           : ')
  cat(pif)
  cat('\n');
  cat('Time steps               : ')
  cat(tif)
  cat('\n');
  cat('Number of records in file:',nRec ,'\n')
  cat('Number of stress periods :',length(pif) ,' \n')
  cat('Number of time steps     :', length(tif),' \n')
#   cat('Stress periods           :', pif ,' \n')
#   cat('Time steps               :', tif,' \n')
  cat('Number of layers         :',NLAY ,' \n')
  cat('Number of Rows           :',NROW ,' \n')
  cat('Number of columns        :',NCOL ,' \n')
  cat('Number of unique labels  :', NLBL,' \n')
  cat('Labels                   :', budlabels,' \n')
}


###############################################################################
##------------------------- READ MODFLOW BINARIES: BUDGET -------------------##
###############################################################################

getValidBudDesc <- function(){
  valid.desc <- c("storage", "constant head", "flow right face",
                    "flow front face", "flow lower face", "wells", "drains",
                    "river leakage")
  return(valid.desc)
}

.Read3dArray <- function(con, nrow, ncol, nl, nbytes) {
  FUN <- function(i) {
    v <- readBin(con, "numeric", n=nrow * ncol, size=nbytes)
    return(matrix(v, nrow=nrow, ncol=ncol, byrow=TRUE))
  }
  return(lapply(seq_len(nl), FUN))
}


bytesOffsetBudget <- function(con){
   valid.desc <- getValidBudDesc()
  for(nskip in 0:4){
  seek(con,where="0",origin="start")  # start position
  if(nskip > 0) skp <- readBin(con, "integer", n=nskip, size=1L)
#   skipBinMarker(con, nskip)
  kstp <- readBin(con, "integer", n=1L, size=4L)
  kper <- readBin(con, "integer", n=1L, size=4L)
  desc <- readBin(readBin(con, "raw", n=16L, size=1L), "character", n=1L)
  desc <- purify(desc)
  if(!desc %in% valid.desc){
      #  cat(paste("variable description not recognized:", desc))
  }else{
#     cat("desc = ",desc, "\n")
#     cat("bytes offset = ",nskip, "\n")
    seek(con,where="0",origin="start")  # start position
    return(as.integer(nskip))
    break
  }
  }
  stop("No bytes offset found\n")
}

firstHdr <- function(con,nskip){    
    skipBinMarker(con, nskip)
#     skp <- readBin(con, "integer", n=1L, size=nskip)
  kstp <- readBin(con, "integer", n=1L, size=4L)
  kper <- readBin(con, "integer", n=1L, size=4L)
  # nbytes <- 4L
  # pertim <- readBin(con, "numeric", n=1L, size=nbytes)
  # totim  <- readBin(con, "numeric", n=1L, size=nbytes)
  desc   <- readBin(readBin(con, "raw", n=16L, size=1L), "character", n=1L)
  desc <- purify(desc)
    nc   <- readBin(con, "integer", n=1L, size=4L)
  nr   <- readBin(con, "integer", n=1L, size=4L)
  nl   <- readBin(con, "integer", n=1L, size=4L)
  skipBinMarker(con, nskip)
#     skp <- readBin(con, "integer", n=1L, size=nskip)
    return(list("kstp"=kstp, "kper" = kper, "desc" = desc, "ncol"=nc, 
                "nrow"= nr, "nlay" = nl))
}

secondHdr <- function(con,nskip,nbytes){
  skipBinMarker(con, nskip)
#     skp <- readBin(con, "integer", n=1L, size=nskip)
    itype  <- readBin(con, "integer", n=1L, size=4L)
  delt   <- readBin(con, "numeric", n=1L, size=nbytes)
  pertim <- readBin(con, "numeric", n=1L, size=nbytes)
  totim  <- readBin(con, "numeric", n=1L, size=nbytes)
  skipBinMarker(con, nskip)
#     skp <- readBin(con, "integer", n=1L, size=nskip)
    return(list("itype"=itype, "delt" = delt, 
                "pertim" = pertim, "totim"=totim))
}

##--- function budget precision f(con, nskip)
budgetPrecision <- function(con, nskip, verbose=FALSE){
  valid.desc <- getValidBudDesc()
  seek(con,where="0",origin="start")  # start position
  hdr1 <- firstHdr(con,nskip)
  compact <- hdr1$nlay < 0

  ncells <- hdr1$ncol * hdr1$nrow * abs(hdr1$nlay)

  fp = ftell(con)    # file position

  precisions <- c('single','double','undefined')
  precisions <- c(4,8,16)
  for(ip in seq_along(precisions)){
  nbytes <- precisions[ip]
  fseek(con,fp,"start")  # reset file pointer
  if(compact){
    ftype <- "compact"
    itype <- secondHdr(con, nskip, nbytes)$itype
    if(itype %in% c(0,1)){
    skipBinMarker(con, nskip)
#     skp <- readBin(con, "integer", n=1L, size=nskip)
    skp <- readBin(con, "numeric", n=ncells, size=nbytes)
#     skp <- readBin(con, "integer", n=1L, size=nskip)
    skipBinMarker(con, nskip)
    }else if(itype %in% c(2)){
    skipBinMarker(con, nskip)
#     skp <- readBin(con, "integer", n=1L, size=nskip)
    nlst <- readBin(con, "integer", n=1L, size=4L)  
    skipBinMarker(con, nskip)
#     skp <- readBin(con, "integer", n=1L, size=nskip)
    if(nlst > 0){
      for(i in 1:nlst){
      skipBinMarker(con, nskip)
#       skp   <- readBin(con, "integer", n=1L, size=nskip)
      icell   <- readBin(con, "integer", n=1L, size=4L)      # ICELL
      val   <- readBin(con, "numeric", n=1L, size=nbytes)  # VAL
      skipBinMarker(con, nskip)
#       skp   <- readBin(con, "integer", n=1L, size=nskip)
      }
    }
    }else{
    stop(' Illegal "itype", must be 0,1 or 2\n')
    }
  }
  hdr1 <- firstHdr(con,nskip)
  desc <- purify(hdr1$desc)
  if(desc %in% valid.desc){
    if(ip == 1){
    if(isTRUE(verbose)){  
      cat('This is a Single Precision Binary Budget File\n')
    }
    pbytes <- 4
    }else if(ip==2){
    if(isTRUE(verbose)){
      cat('This is a Double Precision Binary Budget File\n')
    }
    pbytes <- 8
    }
    if(nskip > 0){
    if(isTRUE(verbose)){
      cat("length of delimiter (bytes):", nskip,"\n")
    }
    }
    return(list("pbytes" = pbytes, "ncol" = hdr1$ncol, 
                "nrow" = hdr1$nrow, "nlay" = hdr1$nlay))
  }
  }
  stop("unable to determine the precision of the budget file\n")
}

budget.info <- function(f){
  con <- file(f, open="rb", encoding="bytes")
  nskip   <- bytesOffsetBudget(con)
  budprec   <- budgetPrecision(con, nskip)
#   ncells   <- budprec$nrow * budprec$ncol * abs(budprec$nlay)
#   nrc    <- budprec$nrow * budprec$ncol
#   pbytes <- budprec$pbytes
  close(con, type="rb")
  return(budprec)
}

getBudget <- function(f){
  valid.desc <- getValidBudDesc()
  con     <- file(f, open="rb", encoding="bytes")
  nskip   <- bytesOffsetBudget(con)
  budprec   <- budgetPrecision(con, nskip)
  ncells   <- budprec$nrow * budprec$ncol * abs(budprec$nlay)
  nrc    <- budprec$nrow * budprec$ncol
  pbytes <- budprec$pbytes
  Nbyt <- flen(con)
  nRecIn     = 0
  nRecOut    = 0
  kperOld    = 0     # period of previous record
  kstpOld    = 0    # kstep  of previous record
  NMAXLBL    = 0    # highest label number on any record
  frewind(con)
  lst <- list()
  scanning <- TRUE
  while(!feof(con)){
  nRecIn <- nRecIn + 1
  hdr1   <- firstHdr(con,nskip)
  kstp   <- hdr1$kstp
  kper   <- hdr1$kper
  desc   <- purify(hdr1$desc)
  nRow   <- hdr1$nrow
  nCol   <- hdr1$ncol
  nl <- hdr1$nlay
  if(feof(con)) cat("file terminer!!!\n")
  if(!desc %in% valid.desc){
    warning(paste("variable description not recognized:", desc,"\n",
      "I will exit now... check if deleting the budget file before 
      running MODFLOW solves the problem!\n"))
  }
  if(kstp > kstpOld){
    cat("\n************ kstp =", kstp, "************\n")
  }
  if(kper > kperOld){
    cat("*----------- kper =", kper, "-----------*\n")
  }
  cat(desc,"\n")
  if(feof(con)){
    break
  }
  if(kper > kperOld || (kper == kperOld && kstp > kstpOld)){
    nRecOut <- nRecOut + 1
    iLb1 <- 1
  }else{
    iLb1 = iLb1 + 1
    if(scanning){
    NMAXLBL = max(NMAXLBL,iLb1)
    }else{
  #     B(nRecOut).label{iLbl} = LABEL;
    }
  }
  if(nl < 0){    # compact
#     cat("compact!\n")
    nl <- abs(nl)
    hdr2   <- secondHdr(con,nskip,pbytes)
    
    itype   <- hdr2$itype
    delt   <- hdr2$delt
    pertim   <- hdr2$pertim
    totim   <- hdr2$totim
    
    if(itype == 5){
    skipBinMarker(con, nskip)
    nval <- readBin(con, "integer", n=1L, size=4L)
    skipBinMarker(con, nskip)
    }else{
    nval <- 1L
    }
    if (nval > 1L) {
    if(nval > 2) stop("nval > 2")
    skipBinMarker(con, nskip)
    ctmp <- readBin(readBin(con, "raw", n=16L, size=1L), 
                    "character", n=nval - 1L)
    ctmp <- purify(ctmp)
    skipBinMarker(con, nskip)
    } else {
    ctmp <- NULL
    } 
    if(itype %in% c(0,1)){
    skipBinMarker(con, nskip)
    d <- .Read3dArray(con, nRow, nCol, nl, pbytes)
    skipBinMarker(con, nskip)
    for (i in seq_along(d)) {
      lst[[length(lst) + 1L]] <- list(d=d[[i]], kstp=kstp, kper=kper,
                      desc=desc, ilay=i, delt=delt,
                      pertim=pertim, totim=totim)
    }
    }else if(itype %in% c(2,5)){
    skipBinMarker(con, nskip)
    nlst   <- readBin(con, "integer", n=1L, size=4L)
    skipBinMarker(con, nskip)
    if(nlst > 0){
      d <- matrix(0, nrow=nlst, ncol=4L + nval)
      colnames(d) <- make.names(c("icell", "layer", "row", "column",
                    desc, ctmp), unique=TRUE)
      for (i in seq_len(nlst)) {
      skipBinMarker(con, nskip)
      d[i, 1] <- readBin(con, "integer", n=1L, size=4L)
      d[i, 4L + seq_len(nval)] <- readBin(con, "numeric", n=nval,
                        size=pbytes)
      skipBinMarker(con, nskip)
      }
      nrc <- nRow * nCol
      d[, "layer"] <- as.integer((d[, "icell"] - 1L) / nrc + 1L)
      d[, "row"] <- as.integer(((d[, "icell"] - (d[, "layer"] - 1L) *
                  nrc) - 1L) / nCol + 1L)
      d[, "column"] <- as.integer(d[, "icell"] - (d[, "layer"] - 1L) *
                    nrc - (d[, "row"] - 1L) * nCol)
      lst[[length(lst) + 1L]] <- list(d=d, kstp=kstp, kper=kper,
                      desc=desc, delt=delt,
                      pertim=pertim, totim=totim)
    }
    }else if(itype %in% 4){
    }else if(itype %in% 3){
    }else{
    stop("unknow itype, must be between 0 and 5")
    }
  }else{    # not compact
    skipBinMarker(con, nskip)
    if(scanning){
      fseek(con, pbytes*ncells);
    }else{
    d <- .Read3dArray(con, nRow, nCol, nl, pbytes)
    }
    skipBinMarker(con, nskip)
  }
  kperOld <- kper
  kstpOld <- kstp
  }
  close(con, type="rb")
  return(lst)
}



###############################################################################
##------------------------- WRITE MODFLOW INPUT FILES ------------------------#
###############################################################################

WriteModflowInputFiles <- function(rs.model, rech = NULL, well = NULL, 
                                   river = NULL, chd = NULL, drain = NULL, 
                                   id = NULL, dir.run = NULL, 
                                   is.convertible = TRUE,
                                   ss.perlen = 0L, tr.stress.periods = NULL, 
                                   ntime.steps = 1L, hnoflo = -999, 
                                   hdry = -888, verbose = TRUE, uni = NULL, 
                                   layConf = NULL, layAvg = NULL, hani = NULL, 
                                   vani = NULL, timeFormat = "%Y%m%d",
                                   header = NULL, wetting=NULL) {
  
  if(is.null(uni)){
    stop("The time and length units must be defined!")
  }else{
    uni[1] <- match.arg(uni[1], c("undefined", "seconds", "minutes", "hours", 
                                "days", "years"))
    uni[2] <- match.arg(uni[2], c("undefined", "feet", "meters", "centimeters"))
    tuni <- which(c("undefined", "seconds", "minutes", "hours", 
                  "days", "years") == uni[1]) -1
    luni <- which(c("undefined", "feet", "meters", "centimeters") == uni[2]) -1
  }
  # lbot <- grep('^lay[[:digit:]]+.bot$',rasterNames) # match layi.bot
  # nl <- length(lbot)    # number of layers
  rasterNames <- names(rs.model)
  nl <- nlay(rs.model)
  
  dir.create(path=dir.run, showWarnings=FALSE, recursive=TRUE)
  if(!is.null(header)){
    header <- paste0("#Groundwater flow model (", Sys.time(), " ",
                     Sys.timezone(), ")")
  }
  iprn <- ifelse(verbose, 3, -1)  # flag for writing to listing file

  # Stress periods
  #  is.transient <- inherits(tr.stress.periods, "Date")
  if(!is.null(tr.stress.periods)){
    is.transient <- TRUE
  }else{
    is.transient <- FALSE
  }

  perlen <- as.integer(ss.perlen)  # stress period length
  nstp   <- 1L    # number of time steps in a stress period
  ss.tr  <- "SS"  # steady-state or transient stress period
  timeID  <- "ss"  # stress period identifier
  if (is.transient) {
    # perlen <- c(0L, as.integer(diff(tr.stress.periods)))
    # nstp <- c(1L, rep(as.integer(ntime.steps), 
    #                   length(tr.stress.periods) - 1L))
    # ss.tr <- c(ss.tr, rep("TR", length(perlen) - 1))
    # timeID <- c(timeID, format(as.Date(head(tr.stress.periods, -1)), 
    #                            timeFormat))
#     perlen <- c(0L, 1L, as.integer(diff(tr.stress.periods)))
    perlen <- c(perlen, 1L, as.integer(diff(tr.stress.periods)))
    nstp <- c(1L, rep(as.integer(ntime.steps), 
                      length(tr.stress.periods) ))
    ss.tr <- c(ss.tr, rep("TR", length(tr.stress.periods) ))
    timeID <- c(timeID, format(as.Date(tr.stress.periods), 
                               timeFormat))
  }

  # Unique unit numbers for output files
  nunit.bud <- 50L  # flow budget
  nunit.hds <- 51L  # hydraulic heads
  nunit.txt <- 52L  # reduced pumping
  
  f <- file.path(dir.run, paste0(id, ".lst"))
  nunit <- 10L
  nam <- data.frame(ftype="LIST", nunit=nunit, fname=basename(f))
  
  #---------------------------------#
  cat("*** Writing MODFLOW Files ***\n")
  cat(length(nstp), "time steps\n")
  cat("units:", uni[1], "and", uni[2], "\n")
  cat("model dimension:", dim(rs.model[[1]])[1], "x", 
      dim(rs.model[[1]])[2], "x", nl, "\n")
  cat("header=", substr(header, 2, nchar(header)) , "\n")
  #---------------------------------#
  
  #==========================
  # Basic file (BAS6)
  cat("Writting BAS6 file...\n")
  f <- file.path(dir.run, paste0(id, ".ba6"))
  nunit <- nunit + 1L
  nam <- rbind(nam, data.frame(ftype="BAS6", nunit=nunit, fname=basename(f)))

  ds.0 <- c(header, "# MODFLOW Basic Package")
  cat(ds.0, file=f, sep="\n", append=FALSE)

  ds.1 <- "PRINTTIME FREE"
  cat(ds.1, file=f, sep="\n", append=TRUE)
  
  # U2DINT
  for (i in 1:nl) {
    ds.2 <- paste("INTERNAL 1 (FREE)", iprn, " # IBOUND layer", i)
    cat(ds.2, file=f, sep="\n", append=TRUE)
    r <- !is.na(rs.model[[paste0("lay", i, ".bot")]])
    r[] <- as.integer(r[])
    m <- format(raster::as.matrix(r), justify = "right", width = 2)
    write.table(m, file=f, append=TRUE, quote=FALSE, row.names=FALSE,
                col.names=FALSE)
  }

  ds.3 <- paste(hnoflo, "HNOFLO")
  cat(ds.3, file=f, sep="\n", append=TRUE)

  for (i in 1:nl) {
    ds.4 <- paste("INTERNAL 1 (FREE)", iprn, " # STRT layer", i)
    cat(ds.4, file=f, sep="\n", append=TRUE)
    m <- raster::as.matrix(rs.model[[paste0("lay", i, ".strt")]])
    m[is.na(m)] <- hnoflo
    write.table(m, file=f, append=TRUE, quote=FALSE, row.names=FALSE,
                col.names=FALSE)
  }

  #==========================
  # Discretization file (DIS)
  cat("Writting DIS file...\n")
  f <- file.path(dir.run, paste0(id, ".dis"))
  nunit <- nunit + 1L
  nam <- rbind(nam, data.frame(ftype="DIS", nunit=nunit, fname=basename(f)))

  ds00 <- c(header, "# MODFLOW Discretization Package")
  cat(ds00, file=f, sep="\n", append=FALSE)
  
  ds01 <- paste(nl, nrow(rs.model), ncol(rs.model), length(perlen), tuni, 
                luni," # NLAY,NROW,NCOL,NPER,ITMUNI,LENUNI")
  cat(ds01, file=f, sep="\n", append=TRUE)
  
  if(is.null(layConf)){
    layConf <- rep(0, nl)
  }else if(length(layConf) != nl){
    stop("The length of 'layConf' (",length(layConf),
         ") is not equal the number of layers (nlay=",nl,")!",sep="")
  }
  ds02 <- paste0(paste(layConf, collapse=" "), " # LAYCBD(NLAY)")
  cat(ds02, file=f, sep="\n", append=TRUE)

  ds03 <- paste("CONSTANT", res(rs.model)[2], " # DELR")
  cat(ds03, file=f, sep="\n", append=TRUE)

  ds04 <- paste("CONSTANT", res(rs.model)[1], " # DELC")
  cat(ds04, file=f, sep="\n", append=TRUE)

  ds05 <- paste("INTERNAL 1 (FREE)", iprn, " # TOP layer 1")
  cat(ds05, file=f, sep="\n", append=TRUE)

  m <- raster::as.matrix(rs.model[["lay1.top"]])
  m[is.na(m)] <- 0.0
  write.table(m, file=f, append=TRUE, quote=FALSE, row.names=FALSE,
              col.names=FALSE)

  for (i in 1:nl) {
    ds06 <- paste("INTERNAL 1 (FREE)", iprn, " # BOTM layer", i)
    cat(ds06, file=f, sep="\n", append=TRUE)
    m <- raster::as.matrix(rs.model[[paste0("lay", i, ".bot")]])
    m[is.na(m)] <- 0.0
    write.table(m, file=f, append=TRUE, quote=FALSE, row.names=FALSE,
                col.names=FALSE)
  }

  ds07 <- data.frame(perlen=perlen, nstp=nstp, tsmult=1, ss.tr=ss.tr)
  write.table(ds07, file=f, append=TRUE, quote=FALSE, row.names=FALSE,
              col.names=FALSE)
              
              
  #===============================
  # Layer-Property Flow file (LPF)
  cat("Writting LPF file...\n")
  f <- file.path(dir.run, paste0(id, ".lpf"))
  nunit <- nunit + 1L
  nam <- rbind(nam, data.frame(ftype="LPF", nunit=nunit, fname=basename(f)))

  ds.0 <- c(header, "# MODFLOW Layer-Property Flow Package")
  cat(ds.0, file=f, sep="\n", append=FALSE)

  if(verbose){
    ilpfcb <- nunit.bud 
  }else{
    ilpfcb <- 0
  }
  # no options as per default in GMS
  ds.1 <- paste(ilpfcb, hdry, 0, " # ILPFCB,HDRY,NPLPF,[Options]")
  #   ds.1 <- paste(ilpfcb, hdry, 0, "CONSTANTCV", "NOVFC",
  #                 "STORAGECOEFFICIENT", " # ILPFCB,HDRY,NPLPF,[Options]")
  cat(ds.1, file=f, sep="\n", append=TRUE)

  if (is.convertible){
    laytyp <- rep(1L, nl)  # convertible = unconfined
  }else{
    laytyp <- rep(0L, nl)  # confined
  }
  ds.2 <- paste0(paste(laytyp, collapse=" "), " # LAYTYP")
  cat(ds.2, file=f, sep="\n", append=TRUE)
  
  if(is.null(layAvg)){
    layAvg <- rep(0L, nl)
  }else if( length(layAvg) != nl){
    stop("The length of 'layAvg' (",length(layAvg),
         ") is not equal the number of layers (nlay=",nl,")!",sep="")
  }
  ds.3 <- paste0(paste(layAvg, collapse=" "), " # LAYAVG")
  cat(ds.3, file=f, sep="\n", append=TRUE)
  if(is.null(hani)){
  # horizontal anisotropy = 1 for all layers (HANI NOT READ)
      hani <- rep(-1L, nl) 
  }else if(length(hani) != nl){
      stop("The length of 'hani' (",length(hani),
       ") is not equal the number of layers (nlay=",nl,")!",sep="")
  }
  ds.4 <- paste0(paste(hani, collapse=" "), " # CHANI")
  cat(ds.4, file=f, sep="\n", append=TRUE)
  
  if(is.null(vani)){
  # VKA is the ratio of horizontal to vertical hydraulic conductivity
  vani <- rep(1, nl)
  # VKA >= 1
  }else  if(length(vani) != nl){
  stop("The length of 'vani' (",length(vani),
       ") is not equal the number of layers (nlay=",nl,")!",sep="")
  }
  ds.5 <- paste0(paste(vani, collapse=" "), " # LAYVKA")
  cat(ds.5, file=f, sep="\n", append=TRUE)
  
  if(is.null(wetting)){
    ds.6 <- paste0(paste(rep(0, nl), collapse=" "), " # LAYWET")
  }else{
    ds.6 <- paste0(paste(rep(1, nl), collapse=" "), " # LAYWET")
  }
  cat(ds.6, file=f, sep="\n", append=TRUE)

  if(!is.null(wetting)){
    ds.7 <- paste0(paste(wetting[1:3], collapse=" "), " # WETFCT IWETIT IHDWET")
    cat(ds.7, file=f, sep="\n", append=TRUE)
  }
  
# no parameters!!!
#   vani <- levels(rs.model[["lay1.zones"]])[[1]]$vani  # TODO(jcf): zones
#   ds.8 <- paste("VERTANISO VANI", vani, "3")
#   cat(ds.8, file=f, sep="\n", append=TRUE)

#   for (i in 1:nl) {
#     ds.9 <- paste(i, "NONE", "ALL")
#     cat(ds.9, file=f, sep="\n", append=TRUE)
#   }

  for (i in 1:nl) {
  # HK (10)
    z <- paste0("lay", i, ".zones")
    
    f.ref <- paste0("hk", i, ".ref")
    fmt <- paste0("OPEN/CLOSE '", f.ref, "' 1.0 '(FREE)' ", iprn, " hk", i)
    cat(fmt, file=f, sep="\n", append=TRUE)
#     cat(ds.13 <- 0, file=f, sep="\n", append=TRUE)
#   cat("*")
    r <- rs.model[[paste0("lay", i, ".hk")]]
  r[is.na(r)] <- hnoflo    # before mv.flag instead of hnoflo
    write.table(raster::as.matrix(r), file=file.path(dir.run, f.ref),
                quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  # HANI (11)
  if(hani[i] <= 0){
    if( paste0("lay", i, ".hani") %in% rasterNames ){
    f.ref <- paste0("hani", i, ".ref")
    fmt <- paste0("OPEN/CLOSE '", f.ref, "' 1.0 '(FREE)' ", iprn, " hani", i)
    cat(fmt, file=f, sep="\n", append=TRUE)
  #     cat(ds.13 <- 0, file=f, sep="\n", append=TRUE)
    r <- rs.model[[paste0("lay", i, ".hani")]]
    r[is.na(r)] <- hnoflo    # before mv.flag instead of hnoflo
    write.table(raster::as.matrix(r), file=file.path(dir.run, f.ref),
          quote=FALSE, row.names=FALSE, col.names=FALSE)
    }else{
      # constant vertical anisotropy = 1
      cat("CONSTANT 1", file=f, sep="\n", append=TRUE)
    }
  }
  
  # VKA (12)
  # if there is a layer for VKA
  if( paste0("lay", i, ".vani") %in% rasterNames ){
    f.ref <- paste0("vani", i, ".ref")
    fmt <- paste0("OPEN/CLOSE '", f.ref, "' 1.0 '(FREE)' ", iprn, " vka", i)
    cat(fmt, file=f, sep="\n", append=TRUE)
#     cat(ds.13 <- 0, file=f, sep="\n", append=TRUE)
    r <- rs.model[[paste0("lay", i, ".vani")]]
    r[is.na(r)] <- hnoflo    # before mv.flag instead of hnoflo
    write.table(raster::as.matrix(r), file=file.path(dir.run, f.ref),
                quote=FALSE, row.names=FALSE, col.names=FALSE)
  }else{
    # constant vertical anisotropy = 1
    cat("CONSTANT 1", file=f, sep="\n", append=TRUE)
  }
  
  if(is.transient){
    f.ref <- paste0("ss", i, ".ref")
    fmt <- paste0("OPEN/CLOSE '", f.ref, "' 1.0 '(FREE)' ", iprn, " ss", i)
    cat(fmt, file=f, sep="\n", append=TRUE)
    r <- deratify(rs.model[[z]], "ss")
    r[is.na(r)] <- hnoflo    # before mv.flag instead of hnoflo
    write.table(raster::as.matrix(r), file=file.path(dir.run, f.ref),
                quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
  #  only if at least one stress period is transient and LAYTYP is not 0.
    if(is.transient && is.convertible){
    f.ref <- paste0("sy", i, ".ref")
    fmt <- paste0("OPEN/CLOSE '", f.ref, "' 1.0 '(FREE)' ", iprn, " sy", i)
    cat(fmt, file=f, sep="\n", append=TRUE)
    r <- deratify(rs.model[[z]], "sy")
    r[is.na(r)] <- hnoflo    # before mv.flag instead of hnoflo
    write.table(raster::as.matrix(r), file=file.path(dir.run, f.ref),
          quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
  
  if(!is.null(wetting) && is.convertible){
#     cat("*")
    cat(paste0("CONSTANT ", wetting[4]), file=f, sep="\n", append=TRUE)
  }
  
  }

  #=================
  # Drain file (DRN)
  if(!is.null(drain)){  
  cat("Writting drain file...\n")
  f <- file.path(dir.run, paste0(id, ".drn"))
  nunit <- nunit + 1L
  nam <- rbind(nam, data.frame(ftype="DRN", nunit=nunit, fname=basename(f)))

  ds.0 <- c(header, "# MODFLOW Drain Package")
  cat(ds.0, file=f, sep="\n", append=FALSE)

  d <- drain[, c("lay", "row", "col", "elev", "cond", "id"), drop=FALSE]

  ds.2 <- paste(nrow(d), nunit.bud, "AUXILIARY id")
  if (!verbose)
    ds.2 <- paste(ds.2, "NOPRINT")
  ds.2 <- paste(ds.2, " # MXACTD,IDRNCB,[Option]")
  cat(ds.2, file=f, sep="\n", append=TRUE)
  ds.5 <- paste(nrow(d), 0, " # ITMP,NP")
  cat(ds.5, file=f, sep="\n", append=TRUE)
  write.table(d, file=f, append=TRUE, quote=FALSE, row.names=FALSE,
        col.names=FALSE)
  if (is.transient) {
    ds.5 <- paste(-1, 0, " # ITMP,NP")
    for (i in 2:length(perlen)) {
    cat(ds.5, file=f, sep="\n", append=TRUE)
    }
  }
  }
  #==================
  # River file (RIV)
  if(!is.null(river)){
  cat("Writting river file...\n")
  f <- file.path(dir.run, paste0(id, ".riv"))
  nunit <- nunit + 1L
  nam <- rbind(nam, data.frame(ftype="RIV", nunit=nunit, fname=basename(f)))

  ds.0 <- c(header, "# MODFLOW River Package")
  cat(ds.0, file=f, sep="\n", append=FALSE)

  n <- nrow(river)
  ds.2 <- paste(n, nunit.bud, "AUXILIARY id")
  if (!verbose)
    ds.2 <- paste(ds.2, "NOPRINT")
  ds.2 <- paste(ds.2, " # MXACTR,IRIVCB,[Option]")
  cat(ds.2, file=f, sep="\n", append=TRUE)

  for (i in seq_along(timeID)) {
    ds.5 <- paste(n, 0L, " # ITMP,NP    STRESS PERIOD", timeID[i])
    cat(ds.5, file=f, sep="\n", append=TRUE)
    d <- river[, c("lay", "row", "col", timeID[i], "cond", "bottom", "id")]
    write.table(d, file=f, append=TRUE, quote=FALSE, row.names=FALSE,
          col.names=FALSE)
  }
  }
  #==========================================
  # Time-variant specified-head package (CHD)
  if(!is.null(chd)){
  cat("Writting specified-head file (CHD)...\n")
  f <- file.path(dir.run, paste0(id, ".chd"))
  nunit <- nunit + 1L
  nam <- rbind(nam, data.frame(ftype="CHD", nunit=nunit, fname=basename(f)))

  ds00 <- c(header, "# MODFLOW River Package")
  cat(ds00, file=f, sep="\n", append=FALSE)
  
  n <- nrow(chd)
  ds02 <- paste(n, nunit.bud, "AUXILIARY abc")
  if (!verbose)
    ds02 <- paste(ds02, "NOPRINT")
  ds02 <- paste(ds02, " # MXACTR,IRIVCB,[Option]")
  cat(ds02, file=f, sep="\n", append=TRUE)

  for (i in seq_along(timeID)) {
    ds05 <- paste(n, 0L, " # ITMP,NP    STRESS PERIOD", timeID[i])
    cat(ds05, file=f, sep="\n", append=TRUE)
    d <- chd[, c("lay", "row", "col", timeID[i], timeID[i], "id")]
#     d <- chd[, c("lay", "row", "col", i, "shead", "ehead", "id")]
    write.table(d, file=f, append=TRUE, quote=FALSE, row.names=FALSE,
          col.names=FALSE)
  }
  }  

  #================
  # Well file (WEL)
  if(!is.null(well)){
  cat("Writting well file...\n")
  f <- file.path(dir.run, paste0(id, ".wel"))
  nunit <- nunit + 1L
  nam <- rbind(nam, data.frame(ftype="WEL", nunit=nunit, fname=basename(f)))

  ds.0 <- c(header, "# MODFLOW Well Package")
  cat(ds.0, file=f, sep="\n", append=FALSE)

  #rech <- as.matrix(cbind(rech[, c("lay", "row", "col", timeID)], id=1))
  #well <- as.matrix(cbind(well[, c("lay", "row", "col", timeID)], id=2))
  #trib <- as.matrix(cbind(trib[, c("lay", "row", "col", timeID)], id=3))
  #misc <- as.matrix(cbind(misc[, c("lay", "row", "col", timeID)], id=4))
  #well <- rbind(rech, well, trib, misc)

  ds.2 <- paste(nrow(well), nunit.bud, "AUXILIARY id")
  if (!verbose)
    ds.2 <- paste(ds.2, "NOPRINT")
  if (is.convertible) {
    ds.2 <- paste(ds.2, "AUTOFLOWREDUCE")
    if (is.transient)
    ds.2 <- paste(ds.2, "IUNITAFR", nunit.txt)
  }
  ds.2 <- paste(ds.2, " # MXACTW,IWELCB,[Option]")
  cat(ds.2, file=f, sep="\n", append=TRUE)

  for (i in seq_along(timeID)) {
    welli <- well[, c("lay", "row", "col", timeID[i], "id"), drop=FALSE]
    #welli <- welli[welli[, timeID[i]] != 0, , drop=FALSE]
    #if (nrow(welli) == 0) next
    ds.5 <- paste(nrow(welli), 0, " # ITMP,NP    STRESS PERIOD", timeID[i])
    cat(ds.5, file=f, sep="\n", append=TRUE)
    write.table(welli, file=f, append=TRUE, quote=FALSE, row.names=FALSE,
          col.names=FALSE)
  }
  }
  
  # RECHARGE file (WEL)
  if(!is.null(rech)){
  cat("NOT writting rech file...\n")
  warnings("package recharge is not yet implemented")
  }
  
  # Sparse Matrix Solver file (SMS)
cat("Writting sparse matrix solver file...\n")
  f <- file.path(dir.run, paste0(id, ".sms"))
  nunit <- nunit + 1L
  nam <- rbind(nam, data.frame(ftype="SMS", nunit=nunit, fname=basename(f)))
  ds.0 <- c(header, "# MODFLOW Sparse Matrix Solver Package")
  cat(ds.0, file=f, sep="\n", append=FALSE)
  iprsms <- if (verbose) 1 else 0
  if (is.convertible) {
    cat("COMPLEX OPTIONS", file=f, sep="\n", append=TRUE)
    ds <- paste(0.001, 0.01, 100, 100, iprsms, 1, 1,
                " # HCLOSE,HICLOSE,MXITER,ITER1,IPRSMS,NONLINMETH,LINMETH")
    cat(ds, file=f, sep="\n", append=TRUE)
  } else {
    ds <- paste(0.001, 0.01, 100, 100, iprsms, 2, 1,
                " # HCLOSE,HICLOSE,MXITER,ITER1,IPRSMS,NONLINMETH,LINMETH")
    cat(ds, file=f, sep="\n", append=TRUE)
    ds <- paste(0.9, 0.0001, 0, 0, 10, 10000, 0.2, 100,
                " # THETA,AKAPPA,GAMMA,AMOMENTUM,NUMTRACK,BTOL,BREDUC,PESLIM")
    cat(ds, file=f, sep="\n", append=TRUE)
    ds <- paste(2, 0, 3, 5, 0, 0, 1, 0.001,
                " # ACL,NORDER,LEVEL,NORTH,IREDSYS,RRCTOL,IDROPTOL,EPSRN")
    cat(ds, file=f, sep="\n", append=TRUE)
  }

  # Output control (OC)
cat("Writting OC file...\n")
  f <- file.path(dir.run, paste0(id, ".oc"))
  nunit <- nunit + 1L
  nam <- rbind(nam, data.frame(ftype="OC", nunit=nunit, fname=basename(f)))

  ds.0 <- c(header, "# MODFLOW Output Control Package")
  cat(ds.0, file=f, sep="\n", append=FALSE)

  ds.1 <- c(paste("HEAD SAVE UNIT", nunit.hds), "COMPACT BUDGET AUXILIARY")
  cat(ds.1, file=f, sep="\n", append=TRUE)

  for (i in seq_along(timeID)) {
    ds.2 <- paste("\nPERIOD", i, "STEP", nstp[i])
    cat(ds.2, file=f, sep="\n", append=TRUE)
    ds.3 <- paste("    ", c("SAVE HEAD", "SAVE BUDGET"))
    cat(ds.3, file=f, sep="\n", append=TRUE)
  }

  # Name file
cat("Writting name file...\n")
  nam <- rbind(nam, data.frame(ftype="DATA(BINARY)", nunit=nunit.bud,
                               fname=paste0(id, ".bud")))
  nam <- rbind(nam, data.frame(ftype="DATA(BINARY)", nunit=nunit.hds,
                               fname=paste0(id, ".hds")))
  if (is.convertible && verbose)
    nam <- rbind(nam, data.frame(ftype="DATA", nunit=nunit.txt,
                                 fname="reduced-pumping.txt"))
  f <- file.path(dir.run, paste0(id, ".nam"))
  write.table(nam, file=f, append=FALSE, quote=FALSE, sep=" ",
              row.names=FALSE, col.names=FALSE)

  invisible(NULL)
}





###############################################################################
##------------------------- WRITE MODPATH INPUT FILES ------------------------#
###############################################################################


setOptionsFlag <- function(ref, new){
  if(!is.null(new)){
  new[names(new) %in% names(ref)]
  ref[names(new)] <- as.integer(new)
  }
  return(ref)
}


writeModpathInputFiles <- function(id=id, dir.run=NULL, optionFlags=NULL,
                                   budgetFaceLabel=NULL,
                  fbud=NULL,  # budget filepath
                  rs.model=NULL, particles = NULL,
                  header=NULL,
                  ReferenceTime=0,
                  Period=1, Step=1, TimeFraction=0,
                  StopTime = 0,
                  GroupCount = 1,  GroupName = "particles",
                  Grid =1, 
                  StartingLocationsFile="jlkjdf.txt",
                  unconfined=TRUE,
                  hnoflo = -999,
                  hdry = -888,
                  verbose = TRUE){   # particle location filename
  
  if(is.null(id) || is.null(dir.run) || is.null(fbud) || is.null(rs.model)){
    stop("'id', 'dir.run', 'fbud', 'rs.model' all must be defined!\n")
  }
  
  iprn <- if (verbose) 3 else -1  # flag for writing to listing file
  
  optionFlagRef <- c("SimulationType" = 2,
          "TrackingDirection" = 1,
          "WeakSinkOption" = 2,
          "WeakSourceOption" = 2,
          "ReferenceTimeOption" = 2,
          "StopOption" = 1,
          "ParticleGenerationOption" = 2,
          "TimePointOption" = 1,
          "BudgetOutputOption" = 1,
          "ZoneArrayOption" = 1,
          "RetardationOption" = 1,
          "AdvectiveObservationsOption" = 1)
  optionFlags <- setOptionsFlag(optionFlagRef,optionFlags)
  
  budgetFaceLabelRef <- c("CONSTANT HEAD" = 6,
           "DRAINS" = 6, 
           "DRAINS (DRT)" = 6,
           "ET SEGMENTS" = 6,
           "ET" = 6,
           "HEAD DEP BOUNDS" = 5,
           "RECHARGE" = 6,
           "RIVER LEAKAGE" = 6,
           "WELLS" = 0)
  budgetFaceLabel <- setOptionsFlag(budgetFaceLabelRef,budgetFaceLabel)
  
  # TO DEFINE!!
  CellBudgetCount <- NULL
  TimePointCount <- NULL
  ReleaseTimeIncrement <- NULL
  TimePoints <- NULL
  TraceFile <- NULL
  TraceID <- NULL
  
  if(!is.null(header)){
  header <- paste0("#Groundwater flow model (", Sys.time(), " ",
          Sys.timezone(), ")")
  }
 
  
  budInfo <- budget.info(fbud)
  NCOL <- budInfo$ncol
  NROW <- budInfo$nrow
  NLAY <- abs(budInfo$nlay)
  
  modflowf0 <- unlist(strsplit(basename(fbud),"[.]"))
  mdflwId <- paste0(modflowf0[1:(length(modflowf0)-1)],collapse=".")

#   optionFlags <- as.integer(optionFlags)
  
  dir.create(path=dir.run, showWarnings=FALSE, recursive=TRUE)
  header <- paste0("# Wood River Valley flow model (", Sys.time(), " ",
                   Sys.timezone(), ")")
     cat("Writting MPSIM file...\n")
     
     ## MPSIM FILE
   f <- file.path(dir.run, paste0(id, ".mpsim"))
   
   #--- COMMENTS
   item0 <- c(header, "# MODPATH simulation file")
   cat(item0, file=f, sep="\n", append=FALSE)
   
   #--- FILENAME AND LISTING FILES
   item1 <- paste0(id, ".mpnam")
   cat(item1, file=f,sep="\n", append=TRUE)
   
   item2 <- paste0(id, ".mplist")
   cat(item2, file=f,sep="\n", append=TRUE)
   
   #--- OPTION FLAGS
   item3 <- paste0(optionFlags,collapse=" ")
   cat(item3, file=f,sep="\n", append=TRUE)
   
   #--- PARTICLE OUTPUT FILENAMES
   item4 <- paste0(id, ".end")
   cat(item4, file=f,sep="\n", append=TRUE)
   
   if(optionFlags["SimulationType"]==2){
    item5 <- paste0(id, ".path")
    cat(item5, file=f,sep="\n", append=TRUE)
   }else if(optionFlags["SimulationType"]==3){
    item6 <- paste0(id, ".ts")
    cat(item6, file=f,sep="\n", append=TRUE)
    if(optionFlags["AdvectiveObservationsOption"]==2){
      item7 <- paste0(id, ".advo")
      cat(item7, file=f,sep="\n", append=TRUE)
    }
   }
   
   #--- REFERENCE TIME
   if(optionFlags["ReferenceTimeOption"]==1){
    if(is.null(ReferenceTime)){
      cat("ReferenceTimeOption set to",0, 
          "but 'ReferenceTime' is not defined!\n")
    }else{
      item8 <- paste0(ReferenceTime)
      cat(item8, file=f,sep="\n", append=TRUE)
    }
   }else if(optionFlags["ReferenceTimeOption"]==2){
    item9 <- paste(Period, Step, TimeFraction, sep=" ")
    cat(item9, file=f,sep="\n", append=TRUE)    
   }
   
   #--- STOPING TIME
   if(optionFlags["StopOption"] == 3){
    item10 <- paste0(StopTime)
    cat(item10, file=f,sep="\n", append=TRUE)    
   }
   
   #--- PARTICLE STARTING LOCATIONS
   # not implemented
   if(optionFlags["ParticleGenerationOption"] == 1){
    stop("ParticleGenerationOption = 1 is not yet implemented!\n")
    item11 <- paste0(GroupCount)
    cat(item11, file=f,sep="\n", append=TRUE)
    
    item12 <- paste0(GroupName)
    cat(item12, file=f,sep="\n", append=TRUE)
    
    Grid <- 1  # must be 1
    item13 <- paste(GroupName)
    cat(item13, file=f,sep="\n", append=TRUE)
    
   }else if(optionFlags["ParticleGenerationOption"] == 2){
    item22 <- paste0(id, ".locations")
    cat(item22, file=f,sep="\n", append=TRUE)
    
    fpart <- file.path(dir.run, item22)
    #--- comments
    pit0 <- c(header, "# MODPATH simulation file")
    cat(pit0, file=fpart, sep="\n", append=FALSE)
    #--- output type
    cat(1L, file=fpart, sep="\n", append=TRUE)
    #--- ground count (how many groups)
    cat(length(particles), file=fpart, sep="\n", append=TRUE)
    for(i in seq_along(particles)){
#         pit2 <- unique(as.character(particles[[i]][,"Label"]))[1]
        pit2 <- paste0("p",i)
        particles[[i]][,"Label"] <- pit2
        #--- group name
        cat(pit2, file=fpart, sep="\n", append=TRUE)
        #--- group particle count
        cat(nrow(particles[[i]]), file=fpart, sep="\n", append=TRUE)
    }
    A0 <- do.call(rbind,particles)
    A <- data.frame("ID"= 1:nrow(A0),"Group"=rep(1:length(particles),
                                                 sapply(particles,nrow)),
            "Grid"=1L, A0)
    write.table(A, file= fpart, append=TRUE,
          quote=FALSE, row.names=FALSE, col.names=FALSE)
    
  
    
  #--- SPECIFIED TIME POINTS
    if(optionFlags["TimePointOption"] == 2 ){
      item23 <- paste0(TimePointCount)
      cat(item23, file=f,sep="\n", append=TRUE)
      item24 <- paste0(ReleaseTimeIncrement)
      cat(item24, file=f,sep="\n", append=TRUE)
    }else if(optionFlags["TimePointOption"] == 3 ){
      item23 <- paste0(TimePointCount)
      cat(item23, file=f,sep="\n", append=TRUE)
      item25 <- paste0(TimePoints)
      cat(item25, file=f,sep="\n", append=TRUE)
    }
   }
   
   #--- NOT IMPLEMENTED: BUDGET OUTPUT
   if( optionFlags["BudgetOutputOption"]==3 ){
    stop("BudgetOutputOption = 3 is not yet implemented!\n")
    item26 <- paste0(CellBudgetCount)
    cat(item26, file=f,sep="\n", append=TRUE)
    for(i in 1:CellBudgetCount){
      item27 <- 0 # GRID LAYER ROW COLUMN of cells for which detailed 
                  # budget information will be generated
      cat(item27, file=f,sep="\n", append=TRUE)
    }
   #--- NOT IMPLEMENTED: BUDGET OUTPUT
   }else  if( optionFlags["BudgetOutputOption"]==4 ){
    stop("BudgetOutputOption = 4 is not yet implemented!\n")
    item28 <- paste0(TraceFile)
    cat(item28, file=f,sep="\n", append=TRUE)
    item29 <- paste0(TraceID)
    cat(item29, file=f,sep="\n", append=TRUE)
   }
   #--- NOT IMPLEMENTED: RETARDATION FACTOR
   if(optionFlags["RetardationOption"] != 1){
    stop("RetardationOption = 1 is not yet implemented!\n")
    item32 <- 0
    cat(item32, file=f,sep="\n", append=TRUE)
    item33 <- 0
    cat(item33, file=f,sep="\n", append=TRUE)
   }
   
   ## MPBAS FILE MODPATH BASIC DATA FILE
   f <- file.path(dir.run, paste0(id, ".mpbas"))
   
   #--- COMMENTS
   item0 <- c(header, "# MODPATH simulation file")
   cat(item0, file=f, sep="\n", append=FALSE)
   
   item1 <- paste(hnoflo, hdry, "HNOFLO HDRY")
   cat(item1, file=f, sep="\n", append=TRUE)
   
   item2 <- paste(length(budgetFaceLabel), "NDefaultIFACE")
   cat(item2, file=f, sep="\n", append=TRUE)
   
   for(i in seq_along(budgetFaceLabel)){
    cat(names(budgetFaceLabel)[i], file=f, sep="\n", append=TRUE)
    cat(paste(budgetFaceLabel[i], "DefaultIFACE"), file=f, sep="\n", 
        append=TRUE)
    
   }
   
   if(unconfined){
    if(NLAY > 1){
      item5 <- paste(c(1,rep(0,NLAY-1)),collapse=" ")
    }else{
      item5 <- "1"
    }
   }else{
    item5 <- paste(rep(0,NLAY),collapse=" ")
   }
   cat(item5, file=f, sep="\n", append=TRUE)
   
   # IBOUND - U2DINT
  for (i in 1:NLAY) {
    item6 <- paste("INTERNAL 1 (FREE)", iprn, " # IBOUND layer", i)
    cat(item6, file=f, sep="\n", append=TRUE)
    r <- !is.na(rs.model[[paste0("lay", i, ".bot")]])
    r[] <- as.integer(r[])
    m <- format(raster::as.matrix(r), justify="right", width=2)
    write.table(m, file=f, append=TRUE, quote=FALSE, row.names=FALSE,
          col.names=FALSE)
  }
  # POROSITY - U2DREL
  for (i in 1:NLAY) {
    f.ref <- paste0("porosity", i, ".ref")
    fmt <- paste0("OPEN/CLOSE '", f.ref, "' 1.0 '(FREE)' ", iprn, 
                  " porosity", i)
    cat(fmt, file=f, sep="\n", append=TRUE)
    r <- rs.model[[paste0("lay", i, ".porosity")]]
    write.table(raster::as.matrix(r), file=file.path(dir.run, f.ref),
          quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  }
   
    ## MPNAM - FILE MODPATH NAME FILE
   f <- file.path(dir.run, paste0(id, ".mpnam"))
   
   item1 <- paste("MPBAS", 70, paste0(id, ".mpbas"))
   cat(item1, file=f, sep="\n", append=FALSE)
   
   item2 <- paste("DIS", 3, paste0(mdflwId, ".dis"))
   cat(item2, file=f, sep="\n", append=TRUE)
   
   item3 <- paste("HEAD", 1, paste0(mdflwId, ".hds"))
   cat(item3, file=f, sep="\n", append=TRUE)
   
   item4 <- paste("BUDGET", 2, paste0(mdflwId, ".bud"))
   cat(item4, file=f, sep="\n", append=TRUE)
} 
   
     

     
################## READ PARTICLE FILES (MODPATH) #################
getSkipLine <- function(con, nc){
  i <- 0
  while(length( newLine <- readLines(con=con,n=1)) > 0){
    if(substr(newLine,1,1) != "#"){
      splLine <- strsplit(gsub("^\\s+|\\s+$", "", newLine), " ")[[1]]
      splLine <- splLine[splLine != ""]
      if(length(splLine) == nc){
      break
      }
    }
#     cat(newLine,"\n")
    i <- i+1
  }
  return(i)
}

readParticles <- function(id, dirRun, ext = NULL, 
                          type = c("loc", "path", "end"), sep = ""){
  type <- match.arg(type,c("loc","path","end") )
  if(type=="loc"){
  fpath <- file.path(dirRun, paste0(id,".locations"))
  fHead <-  c("id","group","grid","lay","row","col",
              "LX","LY","LZ","time","label")  
  }else if(type == "path"){
  fpath <- file.path(dirRun, paste0(id,".path"))
  fHead <- c("id","group","timeIndex","timeStep","time",
             "x","y","z","lay","row",
           "col","grid","lX","lY","lZ","lSegIdx")
  }else if(type == "end"){
  fpath <- file.path(dirRun, paste0(id,".end"))
  fHead <- c("id","group","status","iTime","fTime",
             "iGrid","iLay","iRow","iCol","iFace","iZone",
        "iLX","iLY","iLZ","x0","y0","z0",
             "fGrid","fLay","fRow","fCol","fFace","fZone", 
             "fLX","fLY", "fLZ","x","y","z","label")
  }
  con <- file(fpath,open="r")
  nL <- getSkipLine(con, length(fHead))
  close(con)
#   invisible(readLines(con,i))
#   A <- readLines(con)
  A <- read.table(fpath, skip = nL, header = FALSE, sep = sep)
  names(A) <- fHead
  if(!is.null(ext)){
  if(type == "path" || type == "end"){
    A[,"x"] <- A[,"x"] + ext["xmin"]
    A[,"y"] <- A[,"y"] + ext["ymin"]
#     A[,"z"] <- A[,"z"] + ext["zmin"]
  }
  if(type == "end"){
    A[,"x0"] <- A[,"x0"] + ext["xmin"]
    A[,"y0"] <- A[,"y0"] + ext["ymin"]
#     A[,"z0"] <- A[,"z0"] + ext["zmin"]
  }
  }
  return(A)
}


plotPathXY <- function(P, id = NULL, group = NULL, col = NULL,...){
  
  if(is.null(id) && is.null(group)){
  allid <- unique(P[,"id"])
  allgrp <- unique(P[,"group"])
  if(is.null(col)){
    myCol0 <- plot3D::jet2.col(length(allgrp))
    myCol <- setCol(allid , collim = NULL, col = myCol0)
  }else{
    myCol <- col
  }
#   if(length(col) < length(allid)){
#     col <- col[c(1:length(col),rep(1,length(allid)-length(col)))]
#     cat(length(col))
#   }
  for(i in seq_along(allid)){
    Pid <- P[P[,"id"] == allid[i],]
#invisible(lines(Pid[,c("x")]+xmin(r), Pid[,c("y")]+ymin(r),col=col[i],...))
    invisible(lines(Pid[,c("x")], Pid[,c("y")],col=myCol[i],...))
  }
  
  }else if(!is.null(id)){
  Pid <- P[P[,"id"] == id,]
#   invisible(lines(Pid[,c("x")]+xmin(r), Pid[,c("y")]+ymin(r),...))
  invisible(lines(Pid[,c("x")], Pid[,c("y")],col=col,...))
  }else if(!is.null(group)){
    Pgrp <- P[P[,"group"] == group,]
    allid <- unique(Pgrp[,"id"])
    for(i in seq_along(allid)){
    pp <- Pgrp[Pgrp[,"id"] == allid[i],]
    invisible(lines(pp[,c("x")], pp[,c("y")],col=col,...))
  }
  }
}


# grid <- list(L = 100,   # 0 m to 200 m along x axis
#              W = 40,
#              H = 10,
#              nx = 100,      # number of cells (x axis)
#              ny = 45,      # number of cells (y axis)
#              nz = 20)      # number of cells (z axis)
# prec = 2 # precision for the rounding
# return lay1.top, lay1.bot, ... lay<n>.bot (RasterStack)
modGrid3D <- function(grid, prec=2, fun = NULL,...){
  cat("model with", grid$nx, "(rows) x", grid$ny, "(cols) x",
      grid$nz, "(lays) =", grid$nx * grid$ny * grid$nz, "cells.\n")
  r <- raster(nrows = grid$nx,       ncols = grid$ny, 
              xmn   = grid$W["min"], xmx   = grid$W["max"], 
              ymn   = grid$L["min"], ymx   = grid$L["max"], vals = grid$H)
  names(r) <- "lay1.top"
  if(!is.null(fun)){
    xy <- xyFromCell(r,cellsFromExtent(r, r, expand=FALSE))
    r[] <- fun(x=xy[,1], y=xy[,2],...)
  }
  gwMod <- stack() # initialize a raster stack for model input
  gwMod <- stackRaster(gwMod, r)
  dz <- grid$H/grid$nz
  # bottom layers
  for(i in 1:grid$nz){
    if(i==1){
      r <- gwMod[["lay1.top"]]
    }else{
      r <- gwMod[[paste0("lay",i-1,".bot")]]
    }
    r[] <- round(r[] - dz, prec)
    names(r) <- paste0("lay",i,".bot")
    gwMod <- stackRaster(gwMod, r)
  } 
  return(gwMod)
}

# D = c(L, W, H)  # length, width, heigth
# d = c(dx,dy,dz) # discretization
# prec = 2 # precision for the rounding
# return lay1.top, lay1.bot, ... lay<n>.bot (RasterStack)
modGrid3Dold <- function(D,d,prec=2, fun = NULL,...){
  nRow <- round(D[1]/d[1])
  nCol <- round(D[2]/d[2])
  nLay <- round(D[3]/d[3])
  cat("model with",nRow,"(rows) x", nCol,"(cols) x",nLay,"(lays) =",
      nRow * nCol * nLay,"cells.\n")
  vz <- seq(0, by=d[3], length.out=nLay)
  vx <- seq(0, by=d[1], length.out=nCol)
  vy <- seq(0, by=d[2], length.out=nRow)
  r <- raster(nrows=nRow, ncols=nCol, xmn=min(vx), xmx=max(vx)+d[1], 
              ymn=min(vy), ymx=max(vy)+d[2], vals=max(vz)+d[3])
  names(r) <- "lay1.top"
  if(!is.null(fun)){
    xy <- xyFromCell(r,cellsFromExtent(r, r, expand=FALSE))
    r[] <- fun(x=xy[,1], y=xy[,2],...)
  }
  gwMod <- stack() # initialize a raster stack for model input
  gwMod <- stackRaster(gwMod, r)
  
  # bottom layers
  for(i in 1:nLay){
    if(i==1){
      r <- gwMod[["lay1.top"]]
    }else{
      r <- gwMod[[paste0("lay",i-1,".bot")]]
    }
    r[] <- round(r[]-d[3],prec)
    names(r) <- paste0("lay",i,".bot")
    gwMod <- stackRaster(gwMod, r)
  }  
  return(gwMod)
}

# gwHeads = raster stack 
# fHeads <- file.path(dirRef , paste0("ref3Dss" , ".hds"))
# gwHeads <- get.heads(fHeads, kper = 1:nstp, kstp = 1, r = gwModRef[[1]])
headts <- function(gwHeads, piez, gwMod){

  rnames <- c("lay1.top",paste0("lay", 1:nlay(gwMod),".bot"))
  piez$z <- raster::extract(gwMod[[rnames]], piez$x[,1:2])

  fun1 <- function(x){
    y <- which(x)
    return(tail(y,1))
  }
  fun2 <- function(x){
    y <- which(x)
    return(y[1])
  }
  ltop <- apply(piez$z > piez$x[,3],1, fun1)
  lbot <- apply(piez$z <= piez$x[,3],1, fun2)
  hval <- raster::extract(gwHeads, piez$x[,1:2])
  hts <- list()
  for(i in seq_along(piez$x[,1])){
    ztop <- piez$z[i,ltop[i]] 
    zbot <- piez$z[i,lbot[i]] 
    w <- (ztop - piez$x[i,3])/(ztop - zbot)
    hts[[i]] <- as.vector( 
                  w*hval[i, grep(paste0('^lay',ltop[i]),names(gwHeads))] +
                  (1-w)*hval[i,grep(paste0('^lay',lbot[i]),names(gwHeads))]
                      )
  }
  hts <- do.call("rbind", hts)
  return(hts)
}

# x = layer of hydraulic head
# xmod = rasterStack corresponding to the gw model
saturation <- function(x, xmod, hdry =  -888){
  r <- x > xmod[["lay1.top"]] & x == hdry
  r <- ratify (r)
  rat <- levels(r)[[1]]
  rat$att <- c("Partially Saturated", "Saturated")
  levels(r) <- rat
   grep('^lay[[:digit:]]+.bot$',rasterNames)
  names(r) <- "lay1.saturated"
  return(r)
}

modSection <- function(fhds,kper=1,kstp=1, xmod, vsec=c(1,0), hdry =  -888){
  xhead <- get.heads(fhds,kper=1,kstp=1, r=xmod[[1]], hdry =  -888)
  nl <- nlayers(xhead)
  if(vsec[1]!=0){
    idCells <- cellFromRow(xmod, rownr=vsec[1])
    xap <- xFromCell(xmod, idCells)
  }else{
    idCells <- cellFromCol(xmod, colnr=vsec[2])
    xap <- yFromCell(xmod, idCells)
  }
  z <- matrix(ncol=length(idCells),nrow=nl+1)
  val <- matrix(ncol=length(idCells),nrow=nl)
#   z[1,] <- extract(xmod[[paste0("lay1.top")]], idCells)
  rnames <- c("lay1.top",paste("lay",seq_len(nl),".bot",sep=""))
  z <- raster::extract(xmod[[rnames]], idCells)
  val <- raster::extract(xhead,idCells)
  z <- matrix(z, ncol=length(idCells),nrow=nl+1, byrow=TRUE)
  val <- matrix(val, ncol=length(idCells),nrow=nl, byrow=TRUE)
#   for(i in seq_len(nl)){
#
#     rheads <- xhead[[i]]
#     z[i+1,] <- extract(xmod[[paste0("lay",i,".bot")]], idCells)
#     val[i,] <- extract(xhead[[paste0("lay",i,".head")]], idCells)
#   }
  return(list(x = xap, y = z,z = val))
}

plotPoly <- function(poly, col=rgb(0.2,0.2,0.2,0.4), border=NA,...){
  invisible(sapply(poly, function(X){
    invisible(polygon(X,col=col,border=border,...))
  }))
}

# http://menugget.blogspot.com/2012/04/create-polygons-from-matrix.html#more
r2polygons <- function(r=mat, n=NULL){
  if(missing(r)) stop("Must define matrix 'r'")
  if(missing(n)) stop("Must define at least 1 grid location 'n'")
  poly <- vector(mode="list", length(n))
  xy <- xyFromCell(r,n)
  dx <- res(r)[1]/2
  dy <- res(r)[2]/2
  for(i in seq_along(n)){
    xs <- c(xy[i,1]-dx, xy[i,1]-dx, xy[i,1]+dx, xy[i,1]+dx)
    ys <- c(xy[i,2]-dy, xy[i,2]+dy, xy[i,2]+dy, xy[i,2]-dy)
    poly[[i]] <- data.frame(x=xs, y=ys)
  }
  return(poly)
}

# http://menugget.blogspot.com/2012/04/create-polygons-from-matrix.html#more
m2polygons <- function(r=mat,x, y, n=NULL){
 if(missing(r)) stop("Must define matrix 'r'")
 if(missing(n)) stop("Must define at least 1 grid location 'n'")
 if(missing(x)) x <- seq(0,1,,dim(r)[1])
 if(missing(y)) y <- seq(0,1,,dim(r)[2])
 poly <- vector(mode="list", length(n))
 for(i in seq(n)){
  ROW <- ((n[i]-1) %% dim(r)[1]) +1
  COL <- ((n[i]-1) %/% dim(r)[1]) +1
 
  dist.left <- (x[ROW]-x[ROW-1])/2
  dist.right <- (x[ROW+1]-x[ROW])/2
  if(ROW==1) dist.left <- dist.right
  if(ROW==dim(r)[1]) dist.right <- dist.left
 
  dist.down <- (y[COL]-y[COL-1])/2
  dist.up <- (y[COL+1]-y[COL])/2
  if(COL==1) dist.down <- dist.up
  if(COL==dim(r)[2]) dist.up <- dist.down
 
  xs <- c(x[ROW]-dist.left, x[ROW]-dist.left, 
          x[ROW]+dist.right, x[ROW]+dist.right)
  ys <- c(y[COL]-dist.down, y[COL]+dist.up, 
          y[COL]+dist.up, y[COL]-dist.down)
  poly[[i]] <- data.frame(x=xs, y=ys)
 }
 return(poly)
}



plotSec <- function(xhead, gwMod, vsec=c(0,1),intp=TRUE,
                    col=plot3D::jet.col(101), border=NA,asp=1,...){
#   mySec <- modSection(fhds=fhds,kper=kper,kstp=kstp, gwMod, vsec=vsec)
  # SECTIONING
  nl <- nlayers(xhead)
  if(vsec[1]!=0){
    if(vsec[1] < 1 || vsec[1] > nrow(gwMod)){
      stop("vsec[1] out of range!\n")
    }
    idCells <- cellFromRow(gwMod, rownr=vsec[1])
    xap <- xFromCell(gwMod, idCells)
  }else{
    if(vsec[2] < 1 || vsec[2] > ncol(gwMod)){
      stop("vsec[1] out of range!\n")
    }
    idCells <- cellFromCol(gwMod, colnr=vsec[2])
    xap <- yFromCell(gwMod, idCells)
  }
#   z <- matrix(ncol=length(idCells),nrow=nl+1)
#   val <- matrix(ncol=length(idCells),nrow=nl)
  rnames <- c("lay1.top",paste("lay",seq_len(nl),".bot",sep=""))
  z <- raster::extract(gwMod[[rnames]], idCells)
#   val <- extract(xhead,idCells)
  val <- matrix(raster::extract(xhead,idCells), ncol=length(idCells),
                nrow=nl, byrow=TRUE)
  collim <- range(val)
  if(all.equal(sum(diff(z)),0)){
    # regular grid, just plot the raster
#     cat("lkjL")
    zap <- z[1,]
    zap <- zap[2:length(zap)] - diff(zap)/2
    vval <- t(val)[,nrow(val):1]
    # if there are strange values, they are removed...
    vval[vval < tail(z[1,],1)] <- NA
    image2D(vval,x=xap,y=zap,asp=asp, ...)
    contour2D(vval,x=xap,y=zap, add=TRUE, colkey=FALSE, col="black")
    
  }else{
    z <- matrix(z, ncol=length(idCells),nrow=nl+1, byrow=TRUE)
    pint <- interp(x=rep(xap,each=nl),y=(yy),(val),xo=rev(xap))
    if(intp != TRUE){
      plot(1,1,type="n", xlim=range(xap),ylim=range(z),asp=asp,
           xaxs="i",yaxs="i",bty="n")
      for(k in 1:length(xap)){
        dx <- unique(abs(diff(xap))/2)
        n <- nrow(z)-1
        myCol <- setCol(val[,k], collim=collim, col=col)
        .plotCol(xap[k], z[,k],dx,col=myCol, border=border,...)
        box()
      }
    }else{
      yy <- z[1:(nrow(z)-1),]  + diff(z)/2
      image2D(pint,asp=asp,...)
    }
    contour2D(pint$z,pint$x, pint$y, add=TRUE, colkey=FALSE, col="black")
  }
#   box()
}

modSec <- function(xhead, gwMod, vsec=c(0,1),intp=TRUE){
  nl <- nlayers(xhead)
  if(vsec[1]!=0){
    if(vsec[1] < 1 || vsec[1] > nrow(gwMod)){
      stop("vsec[1] out of range!\n")
    }
    idCells <- cellFromRow(gwMod, rownr=vsec[1])
    xap <- xFromCell(gwMod, idCells)
  }else{
    if(vsec[2] < 1 || vsec[2] > ncol(gwMod)){
      stop("vsec[1] out of range!\n")
    }
    idCells <- cellFromCol(gwMod, colnr = vsec[2])
    xap <- yFromCell(gwMod, idCells)
  }
  rnames <- c("lay1.top", paste0("lay", seq_len(nl), ".bot"))
  z <- raster::extract(gwMod[[rnames]], idCells)
  val <- matrix(raster::extract(xhead, idCells), ncol = length(idCells),
                nrow = nl, byrow = TRUE)
  if( sum(abs(diff(z))) < .Machine$double.eps^0.75){
    # regular grid, just plot the raster
    zap <- as.vector(z[1,])
    zap <- zap[2:length(zap)] - diff(zap)/2
    vval <- t(val)[,nrow(val):1]
    return(list(z = vval, x = xap, y = rev(zap)))    
  }else{
    if(intp != TRUE){
      pp <- vector(mode="list", length(xap))
      #collim <- range(val)
      for(k in 1:length(xap)){
        dx <- unique(abs(diff(xap))/2)
        n <- nrow(z)-1
        # myCol <- setCol(val[,k], collim=collim, col=col)
        pp[[k]] <- list(xleft   = xap[k] - dx, 
                        ybottom = z[2:nrow(z), k], 
                        xright  = xap[k] + dx, 
                        ytop    = z[1:(nrow(z)-1),k],
                        col = val[,k])
#         .plotCol(xap[k], z[,k],dx,col=myCol, border=border,...)
#         box()
      }
      return(pp)
    }else{
      Z <- matrix(z, ncol = length(idCells), nrow = nl+1, byrow = TRUE)
      yy <- Z[1:(nrow(Z)-1),]  + diff(Z)/2
      pint <- interp(x = rep(xap, each = nl), y = yy, z = val, xo = rev(xap),
                     yo = seq(min(yy), to = max(yy), length.out = ncol(yy)))
      return(list(z = pint$z, x = pint$x, y = pint$y))
    }
  }
}

.plotCol <-function(x,y, dx, ...){
rect(x-dx, y[2:length(y)], x+dx, y[1:(length(y)-1)], ...)
}    
    
    
setCol <- function(A , collim = NULL, col = plot3D::jet2.col(n=101)){
  if(is.null(collim)){
    CCY <- (A-min(A,na.rm=TRUE))/(max(A,na.rm=TRUE)-min(A,na.rm=TRUE))
  }else{
    CCY <- (A-collim[1])/(collim[2]-collim[1])
  }
  n <- length(col)
  col[ (CCY)*(n-1) + 1 ] 
}    

dryCells <- function(xhead, xmod){
  nl <- nlay(xmod)
  rnames <- paste("lay",seq_len(nl),".bot",sep="")
  rdry <- stack()
  for(i in 1:nl){
    r <- xhead[["lay1.head"]]
    r[] <- NA
    test <- xhead[[paste0("lay",i,".head")]][] < xmod[[rnames[i]]][] 
    r[test] <- 1
    names(r) <- paste0("lay",i,".dry")
    rdry <- stack(rdry, r)
  }
  return(rdry)
}


gwTable <- function(xhead, r){
#   xhead <- get.heads(fhds,kper=1,kstp=1, r=r[[1]], hdry =  -888)
  return(.gwTable(xhead, r))

}
.gwTable <- function(xhead, xmod){
  r <- xhead[["lay1.head"]]
  rSurfWater <- r
  rSurfWater[] <- NA
  rdry <- r
  rdry[] <- TRUE
  nl <- nlay(xmod)
  rnames <- c("lay1.top",paste("lay",seq_len(nl),".bot",sep=""))
  is.above.land.surface <- xhead[["lay1.head"]][] > xmod[["lay1.top"]][]
  r[is.above.land.surface] <- xmod[["lay1.top"]][is.above.land.surface]
  rSurfWater <- (xhead[["lay1.head"]] - xmod[["lay1.top"]])
  rSurfWater[!is.above.land.surface] <- NA
#   rdry <- rdry & 
  isInLayi <- list()
  for(i in 1:(nl-1)){
    isInLayi <- xhead[[paste0("lay",i,".head")]][] < xmod[[rnames[i]]][] &
          xhead [[paste0("lay",i,".head")]][] > xmod[[rnames[i+1]]][]
    r[isInLayi] <- xhead[[paste0("lay",i,".head")]][isInLayi]
  }
#   i <- i + 1
#   isInLayi <- xhead [[paste0("lay",i,".head")]] < 
#  xmod[[paste0("lay",i-1,".bot")]]
#   r[isInLayi] <- xhead[[paste0("lay",i,".head")]][isInLayi]
  names(r) <- "water.table"
  names(rSurfWater) <- "flooded"
#   namnes(rdry) <- 
  
  return(stack(r, rSurfWater))
}

initialHeads <- function(gwMod, values){
  nLay <- nlay(gwMod)
  r <- gwMod[[1]]
  if(is.numeric(values)){
    for(i in 1:nLay){
      r[] <- values
      r[is.na(gwMod[[paste0("lay",i,".bot")]])] <- NA
      names(r) <- paste0("lay",i,".strt")
      gwMod <- stackRaster(gwMod, r)
    }
    return(gwMod)
  }else if(class(values)[1] == "RasterStack" && nlayers(values) == nLay){
    names(values) <- paste0("lay",1:nLay,".strt")
    for(i in 1:nLay){
      gwMod <- stackRaster(gwMod, values[[i]])
    }
    return(gwMod)
  
  }else if(class(values)[1] == "RasterLayer"){
    for(i in 1:nLay){
      names(values) <- paste0("lay",i,".strt")
      gwMod <- stackRaster(gwMod, values)
    }
    return(gwMod)
    
  }else{
    stop("not yet implemented!\n")
  }
}
    
porosity <- function(gwMod, values){    
  nLay <- nlay(gwMod)
  r <- gwMod[[1]]
  if(is.numeric(values)){
    if(any(values > 1)){
      stop("Porosity values should be between 0 and 1!\n")  
    }
#     p <- stack() # initialize a raster stack for model input
    for(i in 1:nLay){
      r[] <- values
      r[is.na(gwMod[[paste0("lay",i,".bot")]])] <- 0
      names(r) <- paste0("lay",i,".porosity")
      gwMod <- stackRaster(gwMod, r)
    }
    return(gwMod)
  }else{
    stop("not yet implemented!\n")
    # class(values) == raster
  }
}

# co = center
# ro = radius
# n = number of particle on the horizontal circle
# zbot, ztop = vertical limit
# npz = number of "layers" of particles
# total number of particles = n * npz
zylParticles <- function(co=c(0,0),ro=0.5,n=6,zbot,ztop, npz){
  XY <- circPos(co,ro, n)
  zz <- sort(c(zbot,ztop))
#   npz <- 5        # number of particles vertically distributed.
  if(npz > 1){
    pz <- seq(from=zz[1], to=zz[2], by=(zz[2] - zz[1])/(npz-1))
  }else{
    pz <- zz[1] + (zz[2] - zz[1])/(2)
  }
  xyz <- XY[rep(1:nrow(XY),npz),]
  return(cbind(xyz,rep(pz,each=nrow(XY))))
}

circPos <-function(co=c(0,0),ro=0.5,n=6){
  XY <- matrix(nrow=n,ncol=2)
  for(i in seq_len(n)){
    th <- 2*i*pi/n
    rot <- matrix(c(cos(th),sin(th),0,-sin(th),cos(th),0,0,0,1),
                  ncol=3,nrow=3,byrow=TRUE)
    xpos <- rot %*% c(0,1,1)*ro
    XY[i,] <- xpos[1:2] + co
  }
  return(XY)
}

setParticles <- function(gwMod, xyz, pnames = NULL, releaseTime=NULL){
  nLay <- nlay(gwMod)
  r <- gwMod[[1]]
  nxyz <- nrow(xyz)
  if(is.null(pnames)){
    pnames <- paste("p",seq_len(nxyz),sep="")
  }
  idcell <-cellFromXY(r,xyz[,1:2])
  ptclRowCol <- rowColFromCell(r,idcell)
  xyCell <- xyFromCell(r,idcell)
  nptcl1 <- length(idcell)
  px <- round((xyz[,1] %% res(r)[1])/res(r)[1],2)
  py <- round((xyz[,2] %% res(r)[2])/res(r)[2],2)
#   pz <- (xyz[,3] - (xyz[,3] - res(r)[3]/2))/res(r)[3]
  

  if(ncol(xyz) == 2){
    play <- 1
    pz <- 0.5
  }else{
    rnames <- c("lay1.top",paste("lay",seq_len(nLay),".bot",sep=""))
    val <- raster::extract(gwMod[[rnames]],idcell)
    if(any(xyz[,3] > val[,1])){
      stop("The following particles are above the top layer:", 
           which(xyz[,3] > val[,1]),"\n")  
    }
    play <- rep(NA, nxyz)
    pz <- rep(NA, nxyz)
    for(i in seq_len(nxyz)){
#             val <- as.vector(extract(gwMod[[rnames]],idcell[i]))
      vali <- as.vector(val[i,])
      layi <-      which(xyz[i,3] > vali[2:(nLay+1)] & 
                                xyz[i,3] < vali[1:nLay])
            if(length(layi) == 0){
            
            }else{
              play[i] <- layi
              res_z <- vali[play[i]] - vali[play[i]+1]
              pz[i] <- ((xyz[i,3] - vali[play[i]+1]) ) / res_z
      }
    }
  }
  if(is.null(releaseTime)){
    releaseTime <- 0L
  }
  ptcls <- data.frame("Layer" = play,
             "Row" = ptclRowCol[,1],
             "Column" = ptclRowCol[,2],
             "LocalX" = px,
             "LocalY" =  py,
             "LocalZ"  = round(pz,2), 
             "ReleaseTime" = releaseTime,
             "Label" = pnames)
  
  return(ptcls[!is.na(ptcls[,1]),])
}
    
  
    
    
checkModel <- function(xmod, well = NULL, river = NULL, drain = NULL){
  if(!is.null(river)){
    idcell <-cellFromRowCol(xmod[["lay1.bot"]],rownr = river[,"row"], 
                            colnr = river[,"col"])
    lay1bot <- raster::extract(xmod[["lay1.bot"]], idcell)
    if(any(river[,"bottom"] - lay1bot <= 0)){
      cat("*** River bottom below bottom cells for cells (i,j)... ***\n")
      rcfc <- rowColFromCell(xmod[["lay1.bot"]], 
                             idcell[river[,"bottom"] - lay1bot <= 0])
      cat(apply(rcfc,1,paste,sep="\n",collapse=","))  
      cat("\n")
    }
    minH <- apply(river[,7:ncol(river),drop=FALSE],1,min)
    if(any(minH - river[,"bottom"] <= 0)){
      cat("*** River stage below bottom cells for cells (i,j)... ***\n")
      rcfc <- rowColFromCell(xmod[["lay1.bot"]],
                             idcell[minH - river[,"bottom"] <= 0])
      cat(apply(rcfc,1,paste,sep="\n",collapse=","))  
      cat("\n")
    }
  }

}


# CORRECT CHD FRAME
# gwMod = gw model (rasterStack)
# rcCHD = matrix with two columns corresponding to the row and col 
# of the CHD boundaries 
# val = matrix (nrow=number of CHD cells, ncol=number of time step)
# timeID = time index (length = number of time step)
# Check vertical position of the highest CHD value at one position 
# over the time!!!
corCHD <- function(gwMod, rcCHD, val, timeID){
  nl <- nlay(gwMod)
#   rnames <- c("lay1.top",paste("lay",seq_len(nl),".bot",sep=""))
  rnames <- c(paste("lay",seq_len(nl),".bot",sep=""))
  idCells <- cellFromRowCol(gwMod, rcCHD[,"row"],rcCHD[,"col"])
  zl <- as.vector(raster::extract(gwMod[[rnames]],idCells))
#  length(zl)
  Z <- matrix(zl, nrow=nl,byrow=TRUE)
  CHD <- matrix(nrow=nl*length(idCells),ncol=4+length(timeID), 
                dimnames=list(NULL,c("lay","row","col","id",timeID)))
  if(length(dim(val))!=2){
    dim(val) <- c(length(val),1)
  }
  for(k in seq_along(idCells)){
#     lay <- which(max(val[k,]) > Z[,k])
    lay <- which(min(val[k,]) > Z[,k])
  #   CHD[k*nl - nl + 1:nl,"lay"] <- lay
    CHD[nl*(k-1) + lay,"lay"] <- lay
    CHD[nl*(k-1) + lay,"row"] <- rcCHD[k,"row"]
    CHD[nl*(k-1) + lay,"col"] <- rcCHD[k,"col"]
    CHD[nl*(k-1) + lay,c("id")] <- 1L
    CHD[nl*(k-1) + lay,timeID] <- rep(val[k,],each=length(lay))
  }
  CHD <- CHD[!is.na(CHD[,1]),,drop=FALSE]
  return(CHD)
}
    
    
# val <- list(  x = c(18, 14),
#                  y = c(35, 43),
#               ztop = c(7.95,6.9),
#               zbot = c(6.89,3),
#                 id = c(1,2),
#                  q = c(-10000, -5000))
    
setWells <- function(gwMod, val, timeID){
  nl <- nlay(gwMod)
  wellLst <- vector(mode="list",length=length(val$x))
  rnames <- c("lay1.top",paste("lay",seq_len(nl),".bot",sep=""))
  for(k in 1:length(val$x)){
  idCells <- cellFromXY(gwMod,c(val$x[k], val$y[k]))
  zl <- as.vector(raster::extract(gwMod[[rnames]],idCells))
  zl <- round(zl,2)
  aa <- sort(c(val$ztop[k], val$zbot[k]))
  wl <- which(aa[2] >= zl & aa[1] <= zl)
  if(length(wl)==0){
    wl <- which(aa[2] <= zl & aa[1] <= zl)
    AA <- cbind(tail(wl, 1), aa[2]-aa[1])
  }else{
    vv <- c(aa[2], zl[wl], aa[1])
    # layer & proportion
    AA <- cbind(c(wl[1]-1,wl),abs(diff(vv)))
  }
  AA <- AA[AA[,2]!=0,,drop=FALSE]
  if(length(dim(val))==2){
    stop("Not yet implemented for matrix")
  }else{
    qrate <- AA[,2] * val$q[k]/sum(AA[,2])
  }
  if(length(timeID) > 1){
    qrate <- matrix(qrate,nrow=length(qrate),ncol=length(timeID))
  }
  wMat <- cbind(AA[,1],rowFromY(gwMod, val$y[k]), 
                colFromX(gwMod, val$x[k]), val$id[k],qrate)
  wellLst[[k]] <- wMat
  }
  wellFrame <- do.call(rbind, wellLst)
  colnames(wellFrame) <- c(c("lay","row","col", "id"),timeID)
  return(wellFrame)

}

              
# adapted from rasterVis
plotVec <- function(r,scaleSlope=TRUE,dx=1.01,dy=1.01,...){
  SA <- gradh(r)
  sana <- na.omit(SA)
  .arrowsR(sana[,"x"]-0.5*dx*sana[,"dx"], 
           sana[,"y"]-0.5*dy*sana[,"dy"], 
          x1 = sana[,"x"]+0.5*dx*sana[,"dx"], 
           y1 = sana[,"y"]+0.5*dy*sana[,"dy"], ...)
}

gradh <- function(r,scaleSlope=TRUE){
    sa <- terrain(r, opt=c('slope', 'aspect'))
  if (is.logical(scaleSlope) & isTRUE(scaleSlope)){
    sa[["slope"]] <- scale(sa[["slope"]], center = FALSE)
  } else {
    if (is.numeric(scaleSlope)) {
      sa[["slope"]] <- sa[["slope"]]/scaleSlope
    }
  }
  slopex <- sa[["slope"]] * sin(sa[["aspect"]]) 
  slopey <- sa[["slope"]] * cos(sa[["aspect"]])
  SA <- matrix(nrow=ncell(sa),ncol=4)
  SA[,1:2] <- xyFromCell(r, seq_len(ncell(r)))
  SA[,3] <- getValues(slopex)
  SA[,4] <- getValues(slopey)
  colnames(SA) <- c("x","y","dx","dy")
  return(SA)
}

.arrowsR <- function(x0, y0, x1 , y1 ,  length = 0.1, angle = 30,
       code = 2,type= "simple", add = TRUE,...){
    plot3D::arrows2D(x0, y0, x1 , y1,
            length = length, angle = angle,
            code = code,type= type, add = add)
}
              
  
    
