  fun.unitHyd <- function(Q, P = P, nf = 10, mu=0.0001){
    unitHyd(Q, P, nf = nf, mu = mu)
  }
  fun.simHyd <- function(i,uh, P, Q = NULL){
    if(is.null(Q)){
      simHyd(uh[i,], P)
    }else{
      simHyd(uh[i,], P, Q[i,])
    }
  }



splitAt <- function(x, pos) {
  spl <- list()
  k <- 1L
  for(i in  seq_along(pos)){
    spl[[i]] <- x[k:pos[i]]
    k <- pos[i] + 1
  }
  spl[[i+1]] <- x[k:length(x)]
  return(spl)
}

squelette <- function(A, np = 20, fsplit = 3){
  dA <- dist(A)
  pamcl <- cluster::pam(dA, np)
  pc <- pamcl$clustering
  med <- pamcl$id.med
  dMed <- dist(A[med,])
  DD <- as.matrix(dMed)
  X <- A[med, ]
  pts <- numeric(nrow(DD))
  vi <- seq_len(nrow(DD))
  k <- 1
  vk <- c()
  vi <- seq_len(nrow(DD))
  pi <- X[1,,drop = FALSE]
  ksup <- k
  k <- order(DD[k, ])[2]
  pi1 <- X[k, , drop = FALSE]
  vk <- c(vk, ksup)
  for(i in 2:(nrow(DD)-1)){
    pi <- pi1
    ksup <- k
    k <- vi[-vk][order(DD[k, - vk])[2]]
    pi1 <- A[med[k], , drop = FALSE]
    vk <- c(vk, ksup)
  }
  ksup <- k
  vk <- c(vk, ksup)
  
  # scatter2D(A[,1], A[,2], pch = 20)
  # lines(X[vk,], lwd = 2)
  
  # compute distance
  all_dist <- sqrt(apply(diff(A[med,][vk,])^2,1,sum))

  # split
  test <- all_dist > fsplit*sd(all_dist)
  if(sum(test) == 0){
    return(X[vk, ])
  }

  # sort & merge.
  vksplt <- splitAt(vk, which(test))
  vkall <-  vksplt[[1]]
  # 1. first elements - distance with others
  nsp <- (length(vksplt)-1)
  for(i in 1:nsp){
    if(length(vksplt) == 0){
      break
    }
  #   vkall <- c(vkall, vksplt[[i]])
    iht <- c(head(vkall,1), tail(vkall,1))
    vht <- c( sapply(vksplt[-i], head, 1), sapply(vksplt[-i], tail, 1))
    vhtid <- c(1:nsp, - (1:nsp))
    dhi <- apply((t(t(X[vht,]) - X[iht[1],]))^2,1, sum)
    dti <- apply((t(t(X[vht,]) - X[iht[2],]))^2,1, sum)
    if(min(dhi) < min(dti)){
      vkall <- rev(vkall)
      posi <- which.min(dhi)
    }else{
      posi <- which.min(dti)
    }
    id <- vhtid[posi]
    a <- vksplt[-i][[abs(id )]]
    if(id < 0){
      a <- rev(a)
    }
    vkall <- c(vkall, a)
    vksplt <- vksplt[-i][ -abs(id )]
  }
  return(X[vkall,])
}



##---------- FUNCTIONS ----------##

logWrite <- function(fFailure,ok=1,myP){
  cat(paste(c(ok,myP), collapse = "\t"), 
      file = fFailure, sep = "\n", append = TRUE)
}

# # spectral deconvolution with known wavelet
# # convolution model: y = h*x 
# # h and y are known, x is unknown
# # x ~ H^h * Y / (H^h * H + mu)
# deconvFreq <- function(y,h,mu=0.0001){
# 	ny <- length(y)
# 	nh <- length(h)
# 	L  <- ny + ny - 1
# 	H  <- fft(c(h,rep(0,ny-1)))
# 	Y  <- fft(c(y, rep(0,nh-1)))
# 	Re(fft( t(Conj(H))*Y/(t(Conj(H))*H + mu) ,inverse=TRUE))[1:ny]/L
# 	# Re(fft( Y/(H + mu) ,inverse=TRUE))[1:ny]/L
# }

#---------------- CONVOLUTION MODEL --------------------#

forecastModel <- function(prec, rivh, gwh,  weatherFor, convMod){
  timePast <- seq_along(rivh)
  timeFor <- (1L + length(rivh)):length(prec)
  # 0. Unit hydrograph river <-> precipitation
  uh <- unitHyd(y = rivh, x = prec[timePast], nf = convMod$nf, 
                mu = convMod$mu, taper = TRUE, taper2 = TRUE)
  # 1. precipitation forecast
  prec[timeFor] <- simPrec(prec[timeFor], xmin = weatherFor[,1],  
                                          xmax = weatherFor[,2])
  # 2. river stage forecast
  rivSim <- simHyd(uh, P = prec, Q = rivh, f = 1)
  rivSim[timePast] <- rivh
  # 3. Unit hydrograph groundwater <-> river
  uhGW <- t(apply(gwh, 1, fun.unitHyd, nf = convMod$nf, 
                  mu =  convMod$mu, P = rivh ))
  # 4. CHD forecast
  gwSim <- t(sapply(1:nrow(gwh), fun.simHyd, uh = uhGW, 
                      P = rivSim))
  gwSim[,timePast] <- gwh
  return(list(riv = rivSim,
              gw  = gwSim))
}

# spectral deconvolution with known wavelet
# convolution model: y = h*x 
# h and y are known, x is unknown
# x ~ H^h * Y / (H^h * H + mu)
deconvolve <- function(y, h, mu=0.0001, taper = FALSE, nf = 10){
  ny <- length(y)
  nh <- length(h)
  L  <- ny + ny - 1
  H  <- stats::fft(c(h,rep(0,ny-1)))
  Y  <- stats::fft(c(y, rep(0,nh-1)))
  HY <- t(Conj(H))*Y
  HH <- t(Conj(H))*H
  if(isTRUE(taper)){
    nf <- min(nf, L-1)
    tap <- hammingWindow(2*nf -1)
    Rhy <- Re(stats::fft(HY, inverse = TRUE))[,1:nf] * 
              tap[(nf ):(2*nf-1)]
    HY <- fft(c(Rhy, rep(0, L-nf)))
    Rhh <- Re(stats::fft(HH, inverse = TRUE))[,1:nf] * 
              tap[(nf ):(2*nf-1)]
    HH <- fft(c(Rhh, rep(0, L-nf)))
  }
  return( Re(stats::fft( HY/(HH + mu) ,inverse=TRUE))[1:ny]/L)
  # Re(fft( Y/(H + mu) ,inverse=TRUE))[1:ny]/L
}

#' Repeat matrix
#'
#' Repeat a matrix row-wise n times and column-wise m times.
#' 
#' Source
#' A replication of MatLab repmat function!
#' R FOR OCTAVE USERS
#' version 0.4, Copyright (C) 2001 Robin Hankin
#' http://cran.r-project.org/doc/contrib/R-and-octave.txt
#' @name repmat
#' @rdname repmat
#' @export
repmat <- function(A,n,m) {
  kronecker(matrix(1,n,m),A)
}


# pads the edges of an image to minimize edge effects 
# %during convolutions and Fourier transforms. 
# %Inputs %I - image to pad 
# %p - size of padding around image 
# %Output %Ipad - padded image 
# SOURCE: http://matlabgeeks.com/tips-tutorials/how-to-blur-an-image-with-a-
# fourier-transform-in-matlab-part-i/
 # service@matlabgeeks.com i
paddMatrix <- function(I, p1, p2=NULL, zero = FALSE){
  if(is.null(p2)){
    p2 <- p1
  }
  nI <- nrow(I)
  mI <- ncol(I)
  Ipad <- matrix(0, nrow = nI + 2*p1, ncol = mI + 2*p2)
  # middle
  Ipad[(p1+1):(p1+nI),(p2+1):(p2+mI)] <- I
  # top and bottom
  if(zero == FALSE){
    Ipad[1:p1,(p2+1):(p2+mI)] <- repmat(I[1,,drop=FALSE], p1, 1)
    Ipad[(p1+nI+1):(nI+2*p1), (p2+1):(p2+mI)] <- repmat(I[nI,,drop=FALSE], 
                                                        p1, 1)
  }
  if(p2 > 0 && zero == FALSE){
    # left and right
    Ipad[(p1+1):(p1+nI), 1:p2] <- repmat(I[,1,drop=FALSE], 1, p2)
    Ipad[(p1+1):(p1+nI), (p2+mI+1):(mI + 2*p2)] <- 
                            repmat(I[,mI, drop=FALSE],1,p2)
    # corner
    Ipad[1:p1, 1:p2] <- I[1,1]
    Ipad[(p1+nI+1):(nI+2*p1), (p2+mI+1):(mI+2*p2)] <- I[nI,mI]
    Ipad[1:p1,(p2+mI+1):(mI + 2*p2)] <- I[1,mI]
    Ipad[(p1+nI+1):(nI+2*p1), 1:p2] <- I[nI,1]
  }
  return(Ipad)
}


# linear convolution with fft
# a = kernel
convolution <- function(a, b){
	na <- length(a)
	nb <- length(b)
	L <- 2*na + nb - 1
  a0 <- numeric(L)
	b0 <- numeric(L)
	a0[1:na] <- a  #c(a,rep(0,nb-1))
	b0[na + 1:nb] <- b   #c(b, rep(0,na-1))
	y <- Re(fft(fft(a0)*fft(b0),inverse=TRUE))/L
	return(y[na + 1:nb])
}

hammingWindow <- function(L){
	N = L-1
	n <- 0:N
	return(0.54 - 0.46*cos(2*pi*n/N))
}

# estimate unit hydrograph (transfer function)
# from the model Q = U * P
# unitHyd <- function(P,Q, nf = 50, mu=0.01){
	#--- based on the function hydromad::deconvolution.uh.Rd
	# see also Lagged regression models in the frequency domain
	# 4.10 Lagged regression, book R_Time Series Analysis and Its 
	# Applications- with-R-examples
# 	## Calculate PQ and PP cross correlations
# 	nf <- min(nf, length(Q), length(P))
# 	taper <- hammingWindow(2 * nf - 1)
# 	b0 <- ccf(Q, P, lag = nf-1, plot = FALSE, na.action = na.pass)
# 	b <- b0$acf[,1,1] * taper
# 	a0 <- acf(P, lag = nf-1, plot = FALSE, na.action = na.pass)
# 	a <- c(rev(a0$acf[,1,1]), a0$acf[-1,1,1]) * taper
# 	## Deconvolve Q=P*H
# 	taper <- hammingWindow(2 * nf - 1)
# 	uh <- deconvolve(b, a, mu = mu) #* taper[(nf+1):(3*nf-1)]
# 	uh <- uh[1:nf] * taper[(nf):(2*nf-1)]
# model: y = uh * x
unitHyd <- function(y, x, nf = 10, mu=0.0001, taper = TRUE, taper2 = TRUE){
  ny <- length(y)
  nx <- length(x)
  L  <- ny + ny - 1
  H  <- stats::fft(c(x, rep(0, ny - 1)))
  Y  <- stats::fft(c(y, rep(0, nx - 1)))
  HY <- t(Conj(H))*Y
  HH <- t(Conj(H))*H
  if(isTRUE(taper)){
    nf <- min(nf, L-1)
    tap <- hammingWindow(2*nf -1)
    Rhy <- Re(stats::fft(HY, inverse = TRUE))[,1:nf] * 
              tap[(nf ):(2*nf-1)]
    HY <- fft(c(Rhy, rep(0, L-nf)))
    Rhh <- Re(stats::fft(HH, inverse = TRUE))[,1:nf] * 
              tap[(nf ):(2*nf-1)]
    HH <- fft(c(Rhh, rep(0, L-nf)))
  }
  uh <- Re(stats::fft( HY/(HH + mu) ,inverse=TRUE))[1:ny]/L
  if(isTRUE(taper2)){
    tap2 <- hammingWindow(2 * nf - 1)
    uh <- uh[1:nf] * tap2[(nf):(2*nf-1)]
  }
	return(uh)
}

unitHyd2 <- function(P,Q, nf = 10, mu=0.01){
  taper2 <- hammingWindow(2 * nf)
  uh1 <- deconvolve(Q, P, mu = mu)
  uh0 <- uh1[1:nf] * taper2[(nf + 1):(2*nf)]
  return(uh0)
}


# Simulate hydrograph from the convolution
# model Qsim = U * P
# Q is for scaling!!
simHyd <- function(uh, P, Q = NULL, f=0.95){
  Qpred0 <-  convolution(uh, as.numeric(P))
  if(!is.null(Q)){
    Q1 <- as.numeric(Q)
    Qpred0[is.na(Q1)]<-NA
    Qpred0Norm <- Qpred0[seq_along(Q1)]
    Qpred0Norm <- (Qpred0 - min(Qpred0Norm, na.rm=TRUE))/
       (max(Qpred0Norm, na.rm=TRUE) - min(Qpred0Norm, na.rm=TRUE))
    Qpred <- min(Q1, na.rm=TRUE) + 
          Qpred0Norm * f *(max(Q1, na.rm=TRUE)-min(Q1, na.rm=TRUE))
    return(Qpred[seq_along(P)])
  }else{
    return(Qpred0)
  }
}


# precipitation simulation       
simPrec <- function(x, xmin, xmax){
    sapply(seq_along(x), simPrec0, x=x, xmin=xmin, 
                  xmax=xmax)
}
simPrec0 <- function(i, x=0.5, xmin=0, xmax=1, type=c("unif")){
  if(missing(i)){
    i <- 1
  }
#   rsim <- runif(1,xmin[i], xmax[i])
  if(sample(c(0,1),1)==1){
    rsim <- x[i] + abs(rnorm(1, mean=0,sd=abs(xmax[i] - x[i])))
  }else{
    if(x[i] == 0){
      rsim <- 0
    }else{
      rsim <- x[i] - abs(rnorm(1, mean=0,sd=abs(x[i]-xmin[i])))
    }
  }
  rsim <- ifelse(rsim <0, 0, rsim)
  return(rsim)
}


# r = raster identifying position river
# Pend = Modpath output file, end position particles
# Ppath = Modpath output file, path particles
# riverh = riverstage
microbesSim <- function(r = NULL, Pend = NULL,
                        Ppath = NULL, rivh = NULL, 
                        a = 1, b = 0.06, d = 100, span = 0.4,
                        lambda = 0.1){
  nt <- length(rivh)
  timeC <- unique(nt + 1 - Pend[,"iTime"])
  Cw <- numeric(length(timeC))
  # 1. identify particles coming from the river
  test <- extract(r, Pend[,c("x","y")], method = "bilinear")
  vtest <- which(!is.na(test) & !is.na(Pend[,c("z")]))
  if(length(vtest) > 0){
    # 2. compute their path length (3-D)
    Pid <- Pend[vtest, "id"]  # id particles coming from the river
    tRiv <- nt + 1 - Pend[vtest, "fTime"]    # time at river/boundary
    tw <- nt + 1 - Pend[vtest,"iTime"]    # time at well
    dRiv <- numeric(length(vtest))
    for(i in seq_along(Pid)){
      Pi <- Ppath[Pid[i] == Ppath[,"id"],c("x","y","z")]
      Pend[Pend[, "id"] == Pid[i],]
      ni <- seq_len(nrow(Pi) - 1L)
      dRiv[i] <- sum(sqrt(rowSums((Pi[ni, ] - Pi[ni + 1L, ])^2)))
    }
    dRiv <- dRiv[!is.na(dRiv)]
    # 3. microbes concentration in the river for each time step
    Criv <- bacConcRiv2(rivh, a =  a, b = b, d = d, span = span)
    # 4. microbes concentration in the river when each particle reached
    #    the river
    C0 <- signal::interp1(seq_along(rivh), Criv, tRiv,
                          method = c("cubic"), extrap = Criv[1])
    # 5. Mechanical filtration of the microbes by the aquifer
    Cw0 <- 100 * C0 * exp( -0.1 * dRiv)
    nCw <- length(Cw)
    for(k in seq_along(timeC)){
      sel <- tw == timeC[k]
      if( sum(sel) > 0 ){
        Cw[k] <- sum(Cw0[sel])/nCw
      }
    }
  }
  return(Cw)
}
  

# faecal bacteria as function of river stage
bacConcRiv <- function(rivStage, a=1, b=0.1, d=10){
 driv0 <- c(0, -diff(as.numeric(rivStage)))
 driv0[driv0 <=0 ] <- driv0[driv0 <=0 ]/d
 C0 <- 10^(a*driv0 + b)
 return(C0 - min(C0))
}

bacConcRiv2 <- function(rivStage, a=1, b=0.1, d=10, span = 0.5){
 driv0 <- c(0, -diff(as.numeric(rivStage)))
 driv0[driv0 <=0 ] <- driv0[driv0 <=0 ]/d
 C0 <- 10^(a*driv0 + b)
 x <- seq_along(C0)
 C0 <- predict(loess(C0 ~ x, span = span))
 return(C0 - min(C0))
}

# from the model
bacConc <- function(x, xmin,  npt, ptr){
  if(length(dim(x)) == 2){
    Cw <- matrix(ncol=length(ptr),nrow=nrow(x))
    for(k in seq_along(ptr)){
      Cpk <- x[,(k-1)*npt + 1:npt]
      Cw[,k] <- apply(Cpk, 1, estimateConc, xmin)
    }
  }else{
    Cw <- numeric(length(ptr))
    for(k in seq_along(ptr)){
      Cpk <- x[(k-1)*npt + 1:npt]
      Cw[k] <- estimateConc(Cpk, xmin)
    }
  }
  return(Cw)
}




valleyFloor <- function(x,y, a=0.02, b=10.2){
	z <- b - y * a
	return(z)
}

myAsChron <- function(x){
   chron(dates=x, format=c(dates = "d.m.y"))
}

### BUILD MODEL ###
gpHKgwModold <- function(gwMod, meanK_gp, covModelK, sK, cst_mps2mpd ){
  targHK <- xyFromCell(gwMod[["lay1.bot"]], 
                        cellsFromExtent(gwMod[["lay1.bot"]],
                        gwMod[["lay1.bot"]]) )
  KK <- covm(targHK, targHK, covModelK) + diag(sK, nrow(targHK))
  meanlog_gp <- log(meanK_gp) - 0.5 * covModelK$h^2
  paraHydK <-  list("mean"=rep(meanlog_gp, nrow(KK)), "cov"= KK)

  for(i in 1:nlay(gwMod)){
    r <- gwMod[[paste0("lay",i,".bot")]]
    values <- gpSim(paraHydK)
    r[] <- exp(values)*cst_mps2mpd
    r[is.na(gwMod[[paste0("lay",i,".bot")]])] <- NA
    names(r) <- paste0("lay",i,".hk")
    gwMod <- stackRaster(gwMod, r)
  }
  rm(r)
  return(gwMod)
}



gpHKgwMod <- function(gwMod, K_hani, K_hstr, K_vstr, K_nu, K_sd,
                      K_l, K_nug, K_mean, cst_mps2mpd ){
  aniso <- RandomFields::RMangle(angle     = K_hani * pi/180, 
                                lat.angle = 0, 
                                diag      = c(1, K_hstr, K_vstr))
  model <- RandomFields::RMmatern(nu    = K_nu, 
                                  var   = K_sd^2, 
                                  scale = K_l/2, 
                                  Aniso = aniso) +
                                  RandomFields::RMnugget(var = K_nug^2)
  #RFoptions(seed=18051984)
  ya <- yaxis(gwMod)
  vx <- seq(from= ya[1], by = round(mean(diff(ya)),10), 
            length.out = length(ya))
  xa <- xaxis(gwMod)
  vy <- seq(from= xa[1], by = round(mean(diff(xa)),10), 
            length.out = length(xa))
  if(nlay(gwMod) == 1L){
    aniso <- RandomFields::RMangle(angle     = K_hani * pi/180,
                                  diag      = c(1, K_hstr))
    model <- RandomFields::RMmatern(nu    = K_nu, 
                                    var   = K_sd^2, 
                                    scale = K_l/2, 
                                    Aniso = aniso) +
                                    RandomFields::RMnugget(var = K_nug^2)
    simu <- RandomFields::RFsimulate(model, x = vx, y = vy)
  }else{
    #zval <- zaxis(gwMod)
    #zval[1] <- zval[2] - mean(diff(zval[-1]))
    #vz <- seq(from= zval[1], by = round(mean(diff(zval)),10), 
    #        length.out = length(zval))
    zval <- zaxis(gwMod)
    dz <- round(mean(diff(zval)), 10)
    vz <- seq(from = zval[2] - dz/2, 
              by = dz, 
              length.out = nlay(gwMod))
    simu <- RandomFields::RFsimulate(model, x = vx, y = vy, z = vz)
  }
  K0 <- array(as.vector(simu@data)[,1], 
              dim = c(nrow(gwMod), ncol(gwMod), nlay(gwMod)))
  K <-  exp(log(K_mean) + K0)

  for(i in 1:nlay(gwMod)){
    r <- gwMod[[1]]
    r[] <- K[,,i]*cst_mps2mpd
    r[is.na(gwMod[[paste0("lay", i, ".bot")]])] <- NA
    names(r) <- paste0("lay", i, ".hk")
    gwMod <- stackRaster(gwMod, r)
  }
  rm(r)
  return(gwMod)
}


gpHKgwModMultiScale <- function(gwMod, K_hani, K_hstr, K_vstr, K_nu, K_sd,
                      K_l, K_nug, K_mean, cst_mps2mpd, gwModRef ){
  aniso <- RandomFields::RMangle(angle     = K_hani * pi/180, 
                                lat.angle = 0, 
                                diag      = c(1, K_hstr, K_vstr))
  model <- RandomFields::RMmatern(nu    = K_nu, 
                                  var   = K_sd^2, 
                                  scale = K_l/2, 
                                  Aniso = aniso) +
                                  RandomFields::RMnugget(var = K_nug^2)
  #----------------------------------------------
  ya <- yaxis(gwModRef)
  vxref <- seq(from= ya[1], by = round(mean(diff(ya)),10), 
            length.out = length(ya))
  xa <- xaxis(gwModRef)
  vyref <- seq(from= xa[1], by = round(mean(diff(xa)),10), 
            length.out = length(xa))
  if(nlay(gwModRef) == 1L){
    aniso <- RandomFields::RMangle(angle     = K_hani * pi/180,
                                  diag      = c(1, K_hstr))
    model <- RandomFields::RMmatern(nu    = K_nu, 
                                    var   = K_sd^2, 
                                    scale = K_l/2, 
                                    Aniso = aniso) +
                                    RandomFields::RMnugget(var = K_nug^2)
    simuref <- RandomFields::RFsimulate(model, x = vxref, y = vyref)
  }else{
    zval <- zaxis(gwModRef)
    dz <- round(mean(diff(zval)), 10)
    vzref <- seq(from = zval[2] - dz/2, 
                 by = dz, 
                 length.out = nlay(gwModRef))
    simuref <- RandomFields::RFsimulate(model, x = vxref, y = vyref, z = vzref)
  }
   
  K0ref <- array(as.vector(simuref@data)[,1], 
              dim = c(nrow(gwModRef), ncol(gwModRef), nlay(gwModRef)))
  Kref <- exp(log(K_mean) + K0ref)
 
  #----------------------------------------------
  ya <- yaxis(gwMod)
  vx <- seq(from= ya[1], by = round(mean(diff(ya)),10), 
            length.out = length(ya))
  xa <- xaxis(gwMod)
  vy <- seq(from= xa[1], by = round(mean(diff(xa)),10), 
            length.out = length(xa))
  zval <- zaxis(gwMod)
  dz <- round(mean(diff(zval)), 10)
  vz <- seq(from = zval[2] - dz/2, 
            by = dz, 
            length.out = nlay(gwMod))
  
  vxap <- approx(x = vxref, y = seq_along(vxref), xout = vx, rule = 2)$y
  vyap <- approx(x = vyref, y = seq_along(vyref), xout = vy, rule = 2)$y
  vzap <- approx(x = vzref, y = seq_along(vzref), xout = vz, rule = 2)$y
  
  K <- Kref[vxap, vyap, vzap]
 
  # par(mfrow = c(1,2))
  # plot3D::image2D(Kref[,,1], main = "Kref")
  # plot3D::image2D(K[,,1], main = "K-approximate")

  for(i in 1:nlay(gwMod)){
    r <- gwMod[[1]]
    r[] <- K[,,i]*cst_mps2mpd
    r[is.na(gwMod[[paste0("lay", i, ".bot")]])] <- NA
    names(r) <- paste0("lay", i, ".hk")
    gwMod <- stackRaster(gwMod, r)
  }
  rm(r)
  return(gwMod)
}

zonation <- function(gwMod, att.table){
    # delineation of hydrogeologic zones in model layer 1 
  r <- gwMod[["lay1.bot"]]
  r[!is.na(r)] <- 1L
  r <- ratify(r)
  levels(r) <- leftJoin(levels(r)[[1]], att.table,by="ID")
  names(r) <- "lay1.zones"
  gwMod <- stackRaster(gwMod, r)
  for(i in 2:nlay(gwMod)){
    r <- gwMod[["lay1.zones"]]
    names(r) <- paste0("lay",i,".zones")
    gwMod <- stackRaster(gwMod, r)
  }
  return(gwMod)
}

# hrel = relative river stage (between 0 and + something)
# rivH0z = reference elevation of riverstage when hrel = 0
# rivBedz = elevation river bed
# Cr = river bed conductivity
rivGwMod <- function(gwMod, hrel, rivH0z, rivBedz, Cr, timeID){
  cellsRiv <- cellsFromExtent(gwMod[["river"]],trim(gwMod[["river"]])) 
  rivStage <- matrix(hrel, nrow=length(cellsRiv), 
                        ncol=length(timeID), byrow=TRUE)  + rivH0z
  rowColRiv <- rowColFromCell(gwMod[["river"]], cellsRiv)
  riverFrame <- data.frame(lay = 1L,
                          row = rowColRiv[,"row"],
                          col =rowColRiv[,"col"],
                          cond = as.vector(Cr),
                          bottom = as.vector(rivBedz),
                          id=1L)
  riverFrame[,timeID] <- rivStage
  return(riverFrame)
}

### ANALYSIS - DGSA, FPA etc. ###
# find folder names in dirProj with the structure: prfx##### 
# (where ### are numbers)
# return a list with: - id (prfx#####)
#           - nb (#####)
getIdRea <- function(dirProj, prfx = "rea_"){
  idFiles <- list.files(dirProj)
  idRea <- idFiles[grep(paste0("(^", prfx,").([0-9]+$)"), idFiles)]
  idNb <- as.numeric(gsub("[^0-9]", "",idRea))
  return(list(id=idRea, nb=idNb))
}


readHK <- function(dirProj, fName = "HK.rds", ids){
  fHK <- file.path(dirProj, fName)
  n <- length(ids)
  if(!file.exists(fHK)){
    rHK <- as.vector(readRDS(file.path(dirProj, ids[1], fName)))
    HK <- matrix(ncol=length(rHK), nrow=n)
    k <- 0
    while(k < n){
      k <- k + 1
      HK[k,] <- as.vector(readRDS(file.path(dirProj, ids[k], fName)))
    }
    saveRDS(HK, fHK)
  }else{
    HK <- readRDS(fHK)
  }
  return(HK)
}

readRiverFrame <-  function(dirProj, fName = "riverFrame2.rds", ids, 
                            renew = FALSE){
  fRF <- file.path(dirProj, fName)
  n <- length(ids)
  if(!file.exists(fRF) || renew){
    rRF <- (readRDS(file.path(dirProj, ids[1], fName)))
    rRF1 <- as.numeric(rRF[1,7:ncol(rRF)])
    RF <- matrix(ncol=length(rRF1), nrow=n)
    k <- 0
    while(k < n){
      k <- k + 1
      rRF <- (readRDS(file.path(dirProj, ids[k], fName)))
      RF[k,]  <- as.numeric(rRF[1,7:ncol(rRF)])
    }
    saveRDS(RF, fRF)
  }else{
    RF <- readRDS(fRF)
  }
  return(RF)
}


readGWHeads <-  function(dirProj, fName = "gwHeads_lay1.rds", ids, loc,
                            renew = FALSE){
  fHD <- file.path(dirProj, fName)  # file path
  n <- length(ids)
  nloc <- nrow(loc)
  if(!file.exists(fHD) || renew){
    gwHD <- (readRDS(file.path(dirProj, ids[1], fName)))
    piez <- extract(gwHD,loc)
    nt <- ncol(piez)
    bigMat <- matrix(NA, nrow=nloc*n, ncol=nt)
    bigMat[1:nloc,] <- piez
    k <- 1
    while(k < n){
      k <- k + 1
      gwHD <- (readRDS(file.path(dirProj, ids[k], fName)))
      piez <- extract(gwHD,loc)
      bigMat[(k-1)*nloc + 1:nloc,] <- piez
    }
    hds <- splitRowMatToList(bigMat, n=nloc)
    saveRDS(hds, fHD)
  }else{
    hds <- readRDS(fHD)
  }
  return(hds)
}


# matrix formed by alterning 3 rows, split the matrix in three
# x <- matrix(1:18,ncol=2)
# splitRowMatToList(x, n=3)
splitRowMatToList <- function(x, n){
  xspl <- split(t(x), rep(rep(1:n, each = ncol(x)),nrow(x)/n))
  lapply(xspl, matrix, ncol = ncol(x), byrow=TRUE)
}

splitRowMatToList0 <-function(x, n){
  xspl <- split(t(x), rep(rep(1:n, each = ncol(x)),nrow(x)/n))
}

savePrk <- function(dirProj,P, pca1){
  rkHK <- order(pca1)
  frkHK <- file.path(dirProj, "rk_HK.rds")
  saveRDS(rkHK , frkHK)
  
  Paugm <- cbind(P,rkHK)
  fPaugm <- file.path(dirProj, "para_rk_HK.rds")
  saveRDS(Paugm , fPaugm)
  return(Paugm)
  cat(frkHK, "\n")
  cat(fPaugm, "\n")
}

estimateConc <- function(x, xmin){
  nx <- length(x)
  x <- x[!is.na(x) & x >= xmin] # concentration above detection threshold
#   return( (length(x)/nx) * sum(x) )
  return( sum(x)/nx )
}



#' Interpolate tracer concentration to get nice curves
#'
#' Use a Gaussian smoothing
#' @para w length-one numeric vector defining the window length of the smoothing
#' @para n length-one integer vector defining the number of monitoring locations
#'          (m must be a multiple of n)
#' @para D (n x m) data matrix (n = # of realisation, m = # of sample)
tracerInterp <- function(D, n, w = 5){
  if(is.null(dim(D))){
    dim(D) <- c(1,length(D))
  }
  ntD <- ncol(D)/n
  Dnew <- list()
  for(i in 1:n){
    sel <- (i-1)*ntD + 1:ntD
    D0 <- D[, sel, drop=FALSE]
    D00 <- cbind(matrix(0,nrow=nrow(D0), ncol=50), D0)
    Dnew[[i]] <- t(apply(D00, 1, mmand::gaussianSmooth, 5))[,51:ncol(D00)]
  }
  return(Dnew)
}

#' Split the observation into a list
splitObs <- function(obs, n){
  if(is.null(dim(obs))){
    ntD <- length(obs)/n
    obsNew <- list()
    for(i in 1:n){
      sel <- (i-1)*ntD + 1:ntD
      obsNew[[i]] <- obs[sel]
    }
  }else{
    ntD <- ncol(obs)/n
    obsNew <- list()
    for(i in 1:n){
      sel <- (i-1)*ntD + 1:ntD
      obsNew[[i]] <- obs[, sel]
    }
  }
  return(obsNew)
}

cumsump <- function(x, name = NULL){
  if(is.null(name)){
    cumsum(x)/sum(x)
  }else{
    cumsum(x[[name]])/sum(x[[name]])
  }
}

# robust PCA wrapper for pcaPP::PCAproj
robpca <- function(x, k = NULL, scale = NULL,  maxit = 5, 
                  center = l1median_NLM, method = c("mad", "sd","qn") ){
  if(is.null(k)){
    k <- ncol(x)
  }
  P <- PCAproj(x, k = k, scale = scale, maxit = maxit, center = center,
                method = method)
  P$rotation <- P$loadings[]
  P$loadings <- NULL
  P$x <- P$scores[]
  P$scores <- NULL
  if(is.null(scale)){
    P$scale <- FALSE
  }
  return(P)
}

# robust CCA wrapper for pcaPP::PCAproj
robcca <- function(X, Y, k = NULL, method = c("spearman"), standardize = FALSE,
                    useL1Median = FALSE){
  if(is.null(k)){
    k <- min(ncol(X), ncol(Y))
  }
  P <- CCAproj( X, Y, k = k, method = method, standardize = standardize, 
          useL1Median = useL1Median)
  P$xcoef <- P$A
  P$ycoef <- P$B
  P$scores <- comput(X,Y, P)
  return(P)
}

##################### PREDICTION-FOCUSED ANALYSIS ######################
# # H = simulated forecast
# # D = simulated data
# # Dobs = observed data
# # nFor = number of sample of the posterior p(H | D, Dob)
# pfa <- function(H, D, Dobs, thr = 99, nFor = 1000){
#   # forecast
#   H_PCA <-  prcomp(H, center = TRUE, scale = FALSE)
#   HcumVar <- cumsum(100*(H_PCA$sdev)^2/sum((H_PCA$sdev)^2))
#   # data
#   D_PCA <-  prcomp(D, center = TRUE, scale = FALSE)
#   DcumVar <- cumsum(100*(D_PCA$sdev)^2/sum((D_PCA$sdev)^2))
#   # observed data
#   Dobs_PCA <- projPCA(x = Dobs, PCA = D_PCA)
#   nD <- tail(which(DcumVar < thr),1)
#   nH <- tail(which(HcumVar < thr),1)
#   if(length(nH) == 0){
#     nH <- 1
#   }
#   X <- D_PCA$x[,seq_len(nD), drop = FALSE]
#   Y <- H_PCA$x[,seq_len(nH), drop = FALSE]
#   if(nH == 1){
#     post <- blrpred(x = X, y = Y, xobs = Dobs_PCA[,1:nD,drop=FALSE],
#                     N = nFor)
#     forH_PCA <- post
#     Dobs_CCA <- NULL
#     forH_CCA <- NULL
#     CCA <- list(X=X, Y = Y, Xobs =  Dobs_PCA, Ypred = forH_PCA, n = nD)
#   }else{
#     CCA <- cc(X, Y)
#     Dobs_CCA <- projCCA(x = Dobs_PCA[,1:nD,drop=FALSE], 
#                         CCA, center = colMeans(X))
#     #--- Posterior mean and covariance
#     # post <- postGaussian(X = CCA$scores[[1]], Y = CCA$scores[[2]], 
#     #                       xobs = Dobs_CCA, outliers = TRUE)
#     # forH_CCA <- MASS::mvrnorm(n = nFor, mu = post$mean, Sigma = post$cov)
#     forH_CCA <- matrix(nrow = nFor, ncol = length(CCA$cor))
#     for(k in seq_along(CCA$cor)){
#       forH_CCA[,k] <- blrpred(y = CCA$scores[[2]][,k], 
#                       x = CCA$scores[[1]][,k], 
#                       xobs = Dobs_CCA[,k], N = nFor)
#       # forH_CCA[,k] <- mvtnorm::rmvt(n = nFor, delta = post$delta, 
#       #                              sigma = post$sigma, df = post$df)
#     }
#     # post <- lmpred(X = CCA$scores[[1]], Y = CCA$scores[[2]], 
#     #                    xobs = Dobs_CCA, npred = nFor, outliers = TRUE)
#     #forH_CCA <- post$pred
#     #---- back projection
#     forH_PCA <- projCCA(forH_CCA, CCA, type="y", inverse =TRUE, 
#                       center =  colMeans(Y))
#   }
#   #----- 
#   forH <- projPCA(forH_PCA, H_PCA, inverse = TRUE)
#   return(list(CCA = CCA, obsCCA = Dobs_CCA, predCCA = forH_CCA, pred = forH))
# }

# @param Htransf either NULL, 'log' or 'nScore'
PFA <- function(D1, D2, H,D1obs, D2obs, cumsd = 99, cumsd0=100,
                what = "heads", nFor =1000, 
                transf = FALSE, nD = NULL, nH = NULL, test  = FALSE, 
                nperm = NULL, Htransf = NULL, HtransfPar = 0.0001){
  H0 <- H
  ## FORECAST
  if(!is.null(Htransf) && Htransf == "log"){
    H <- log(H + HtransfPar[1])
  }else if(!is.null(Htransf) && Htransf == "sinh"){
    # Inverse hyperbolic sine (IHS) transformation
    H <- log(H*HtransfPar[1] + sqrt((H * HtransfPar[1])^2 + 1))/HtransfPar[1]
  }else if(!is.null(Htransf) && Htransf == "nScore"){
    H_nscores <- apply( H, 2, nScoreTrans)
    H <- do.call(cbind, lapply(H_nscores,function(xl) xl$x) )
    Hns_tbl <- lapply(H_nscores,function(xl) xl$table)
  }
  H_PCA <-  prcomp(H,center=TRUE, scale=FALSE)
  ## DATA
  if(what=="tracer" || what == "both"){
    D1_fpca <- multiFPCA(D1, nH = 100, nbasis = 100, norder = 3, 
                        centerfns = TRUE,
                        scaling = "global")
    nHarmD1 <- max(sapply(lapply(lapply(D1_fpca, cumsump, name="values"), '<=', 
                  cumsd0/100), sum))
    nHarmD1 <- min(nHarmD1, length(D1_fpca[[1]]$varprop))
    #--- Project observations on fPCA
    D1obs_fpca <- projfPCA(x = D1obs, fPCA = D1_fpca, nH = nHarmD1)
    #--- 4. PCA for all wells
    D1_fpcaScores <- fPCAScores(D1_fpca, nH = nHarmD1)
  }
  if(what=="heads" || what == "both"){
    D2_fpca <- multiFPCA(D2, nH = 30, nbasis = 30, norder = 3, 
                          centerfns = TRUE,scaling = "global")
    nHarmD2 <- max(sapply(lapply(lapply(D2_fpca, cumsump, name="values"), '<=', 
                  cumsd0/100), sum))
    nHarmD2 <- min(nHarmD2, length(D2_fpca[[1]]$varprop))
    nHarmD2 <- min( length(D2_fpca[[1]]$varprop))
    cat(nHarmD2)
    #--- Project observations on fPCA
    D2obs_fpca <- projfPCA(x = D2obs, fPCA = D2_fpca, nH = nHarmD2)
    #--- 4. PCA for all wells
    D2_fpcaScores <- fPCAScores(D2_fpca, nH = nHarmD2)
  }
  if(what == "both"){
    D_f <- cbind(D1_fpcaScores, D2_fpcaScores)
    D_f <- cbind(D1_fpcaScores, D2_fpcaScores)
  }else if(what == "heads"){
    D_f <- D2_fpcaScores
  }else if(what == "tracer"){
    D_f <- D1_fpcaScores
  }

  D_fPCA <- prcomp(D_f, center = FALSE, scale = FALSE)
  # OBS > PCA on fPCA
  if(what == "both"){
    Dobs_f <- c(do.call(c, D1obs_fpca), do.call(c, D2obs_fpca))
  }else if(what == "heads"){
    if(is.null(dim(D2obs_fpca[[1]]))){
      Dobs_f <- do.call(c, D2obs_fpca)
    }else{
      Dobs_f <- do.call(cbind, D2obs_fpca)
    }
  }else if(what == "tracer"){
    if(is.null(dim(D1obs_fpca[[1]]))){
      Dobs_f <- do.call(c, D1obs_fpca)
    }else{
      Dobs_f <- do.call(cbind, D1obs_fpca)
    }
#     Dobs_f <- do.call(c, D1obs_fpca)
  }

  Dobs_fPCA <- projPCA(x = Dobs_f, PCA = D_fPCA)
  
  #--- normal quantile transformation
  X <- D_fPCA$x
  Y <- H_PCA$x
  if(isTRUE(transf)){
    X_nscores <- apply( X, 2, nScoreTrans)
    Xt <- do.call(cbind, lapply(X_nscores,function(xl) xl$x) )
    Xtbl <- lapply(X_nscores,function(xl) xl$table)
    FUN <- function(i, x, tbl){
      nScoreTrans(x[,i], tbl = tbl[[i]])$x
    }
    Dobs_c <- sapply(1:length(D_fPCA$sdev), FUN, Dobs_fPCA, tbl = Xtbl)
    if(is.null(dim(Dobs_c))){
      dim(Dobs_c) <- c(1,length(Dobs_c))
    }
    Y_nscores <- apply( Y, 2, nScoreTrans)
    Yt <- do.call(cbind, lapply(Y_nscores,function(xl) xl$x) )
    Ytbl <- lapply(Y_nscores,function(xl) xl$table)
#     Yt <- Y
  }else{
    Dobs_c <- Dobs_fPCA
    Xt <- X
    Yt <- Y
  }
  
  ## CCA
  DcumVar <- cumsum(100*(D_fPCA$sdev)^2/sum((D_fPCA$sdev)^2))
  HcumVar <- cumsum(100*(H_PCA$sdev)^2/sum((H_PCA$sdev)^2))
  nX <- which(DcumVar >= cumsd)[1]
  nX <- ifelse(is.na(nX), length(DcumVar), nX)
  nY <- which(HcumVar >= cumsd)[1]
  nY <- ifelse(is.na(nY), length(HcumVar), nY)
  if(is.null(nD)){
    nD <- nX
  }else{
    nD <- min(length(D_fPCA$sdev),nD)
  }
  if(is.null(nH)){
    nH <- nY
  }else{
    nH <- min(ncol(H),nH)
  }
  if(test == TRUE){
    nMaxD <- length(D_fPCA$sdev)
    nMaxH <- ncol(H)
    myRes <- numeric(max(nMaxD,nMaxH)-2)
    for(i in 3:max(nMaxD, nMaxH)){
      X <- Xt[,2:min(i, nMaxD), drop = FALSE]   # D
      Y <- Yt[,2:min(i, nMaxH), drop = FALSE]
      CCA <- cc(X, Y)
      Yb <- CCA$scores[[2]] %*% pinv(CCA$ycoef)
      Cwb <- projPCA(Yb, H_PCA, inverse = TRUE, comp = 1:ncol(Yb))
      res <- (Cwb- H)
      myRes[i-2] <- mean(abs(res)) 
    }
    test <- 100*abs(diff(myRes))/(max(myRes)-min(myRes))
     i <- max(which(test > 1)) + 2
     nD <- min(i, nMaxD)
     nH <- min(i, nMaxH)
  }
  X <- Xt[,1:nD, drop = FALSE]   # D
  Y <- Yt[,1:nH, drop = FALSE]   # H

  ### USE vegan::CCorA() with nperm=1000!!!! (cf. R_Regression avec R)
  ### CHECK: Noise and Outlier Filtering in Heterogeneous Medical Data Sources
  ### CHEKC: RGCCA::rgcca
  ### CHEKC: ade4::coinertia
  
  
  CCA <- cc(X, Y)

  if(!is.null(nperm)){
    pval <- permCCA(X,Y, nperm=100)
    pasym <- p.asym(CCA$cor,nrow(X), ncol(X), ncol(Y))
  }else{
    pval <- NULL
    pasym <- NULL
  }
  # OBS > CCA
  Dobs_CCA <- projCCA(x = Dobs_c[,1:nD,drop=FALSE], CCA, center = colMeans(X))

  #--- Posterior mean and covariance
  post <- postGaussian(X = CCA$scores[[1]], Y = CCA$scores[[2]], 
                      xobs = Dobs_CCA, outliers = TRUE)
  FUN <- function(i, n, mu, Sigma){
    mvrnorm(n=n, mu=mu[i,], Sigma = Sigma)
  }
  forCCA <- lapply(1:nrow(post$mean), FUN , n = nFor, mu=post$mean,
          Sigma= post$cov)              
#   forCCA <- mvrnorm(n = nFor, post$mean, post$cov)

  #### Statistical Methods for Multivariate Outlier Detection
  #### --> use mahalanobis distance
  forPCA  <- lapply(forCCA, projCCA, CCA, type="y", inverse =TRUE, 
                    center =  colMeans(Yt))
#   forPCA <- projCCA(forCCA, CCA, type="y", inverse =TRUE, 
#                     center =  colMeans(Yt))
  #--- backtransform
   # residuals reconstruction
  Yb <- CCA$scores[[2]] %*% pinv(CCA$ycoef)
  
  if(isTRUE(transf)){
    FUN1 <- function(i, x, tbl){
      nScoreTrans(x[,i], inverse = TRUE, tbl = tbl[[i]])
    }
    Yb <- sapply(1:ncol(Yb), FUN1, Yb,  tbl = Ytbl)
    FUN2 <- function(x, tbl = Ytbl){
      sapply(1:ncol(x), FUN1, x, tbl = Ytbl)
    }
    forPCAb <- lapply(forPCA, FUN2, tbl = Ytbl)
#     forPCAb <- sapply(1:ncol(forPCA), FUN1, forPCA, tbl = Ytbl)
#     forPCAb <- forPCA
  }else{
    forPCAb <- forPCA
  }
  
  # residuals reconstruction
  Cwb <- projPCA(Yb, H_PCA, inverse = TRUE, comp = 1:ncol(Yb))
  
  forCw <- lapply(forPCAb, projPCA, H_PCA, inverse = TRUE)
#   forCw <- projPCA(forPCAb, H_PCA, inverse = TRUE)
  
  if(!is.null(Htransf) && Htransf == "log"){
    forCw <- lapply(forCw, exp)
    forCw <- lapply(forCw, '-', HtransfPar[1])
#     forCw <- exp(forCw) - HtransfPar[1]
    Cwb <- exp(Cwb) - HtransfPar[1]
  }else if(!is.null(Htransf) && Htransf == "sinh"){
    forCw <- lapply(forCw, '*' , HtransfPar[1])
    forCw <- lapply(forCw, sinh)
    forCw <- lapply(forCw, '/', HtransfPar[1])
     Cwb <- sinh(Cwb*HtransfPar[1]) / HtransfPar[1]
  }else if(!is.null(Htransf) && Htransf == "nScore"){
#     forCw <- sapply(1:ncol(forCw), FUN2, forCw,  tbl = Hns_tbl)
    forCw <- lapply(forCw, FUN2, tbl = Hns_tbl)
    Cwb <- sapply(1:ncol(Cwb), FUN1, Cwb, tbl = Hns_tbl)
  }
  if(length(forCw) == 1){
    forCw <-  forCw[[1]]
    forCCA <- forCCA[[1]]
  }
  
    res <- (Cwb- H0)
  resReconstruction <- summary(abs(res)) 
  
  return(list(pred = forCw, predCCA = forCCA, CCA = CCA, 
              obsCCA = Dobs_CCA, outliers = post$outliers, nD = nD, 
              nH = nH, nDvar = nX, nHvar = nY, pval = pval, pasym = pasym,
              res = resReconstruction, backTrans = Cwb))
}


margLik <- function(X, x, kernel = "normal", canonical = FALSE, 
                    gridsize = 401L, truncate = TRUE){
  ksPred <- bkde(X, kernel = kernel, canonical = canonical, 
                 gridsize = gridsize, truncate = truncate)
  df <- approxfun(ksPred, rule=2)
  return(df(x))
}

# self-consistency
# marginal likelihood
checkPFA <- function(iObs = NULL, D1, D2, H,  ...){
   vn <- 1:nrow(H)
  if(is.null(iObs)){
    iObs <- vn
  }
#   vn <- vn[-iObs]
  D1 <- lapply(D1,function(x, n) x[n,], vn)
  D1obs <- lapply(D1,function(x, n) x[n,], iObs)
  D2 <- lapply(D2,function(x, n) x[n,], vn)
  D2obs <- lapply(D2,function(x, n) x[n,], iObs)
#   H <- H0[vn,]
  Hobs <- H[iObs,]
  PPAll <- PFA(D1, D2, H, D1obs, D2obs, ...)
  # posterior prediction
  postPred_tracer <- do.call(rbind,PPAll$pred)
  FUN <- function(i, x, a, b = 9){
    ML0 <- 0
    for(k in 1:ncol(x[[i]])){
        mltest <- margLik(x[[i]][,k], a[i,k])
        mltest <- ifelse(mltest < 10^(-300), 10^(-300), mltest)
       ML0 <- ML0 + log( mltest)
    }
#     margLik(x[[i]][,b], a[i,b])
    return(ML0)
  }
  # prediction performance for each observations
  MLpred <- sapply(seq_along(PPAll$pred), FUN, PPAll$pred, Hobs)
  return(list(pred=postPred_tracer, logML = MLpred, nD = PPAll$nD, 
              nH = PPAll$nH))
}

# 
# plotREG <- function(PP, colvar = NULL, col = NULL, pch = 1, colvar2 = NULL,
#                     col2="firebrick3", pch2 = 1, xlab = "data", 
#                     ylab = "forecast"){
# #    if(is.null(colvar) && !is.null(PP$outliers)){
# #     colvar <- !PP$outliers
# #    }
#   par.default <- par(no.readonly=TRUE)
#    if(is.null(colvar2)){
#     colvar2 <- NULL
#    }
#    par(mfrow=c(2,2), cex =1, cex.lab=1, mai=c(1.2,1.2,1.2,0.4))
#     ni <- min(4, (PP$CCA$n))
#     for(i in 1:ni){
#       scatter2D(PP$CCA$X[,i], PP$CCA$Y[,1], 
#                 colvar = colvar, pch = pch,col = col,
#                 xlim = range(c(PP$CCA$X[,i], PP$CCA$Xobs[1,i])),
#                 ylim = range(c(PP$CCA$Y[,1], PP$CCA$Ypred[,1])),
#                 xlab=xlab, ylab=ylab, asp=NA, colkey=FALSE)
#               #  main=paste0(i," canonical variables\n(correlation = ", 
#               #  round(PP$CCA$cor[i],2),")"))
#       abline(v=PP$CCA$Xobs[1,i], col="chartreuse3", lwd=4)
#       scatter2D(rep(PP$CCA$Xobs[1,i], length( PP$CCA$Ypred[,1])), 
#                   PP$CCA$Ypred[,1], 
#                 pch=pch2, col=col2, cex=1.5, add = TRUE, colvar = colvar2,
#                 colkey = TRUE)
#       legend("topleft", c("prior","posterior","obs. data"), 
#             col=c("dodgerblue2","firebrick3","chartreuse3"), 
#           pch=c(pch, pch2, NA), lwd = c(NA,NA,2), bty="n", cex=1)
#             grid()
#     }
# 	par(par.default)
#     # dev.off()
# }


# CairoPNG(filename = file.path(paste0(dirProj, "_00_", what, 
#           "_CCA_with_samples.png")), 
#          width = 900, height = 900, pointsize = 18)
# plotCCA <- function(PP, colvar = NULL, col = NULL, pch = 1, colvar2 = NULL,
#                     col2="firebrick3", pch2 = 1){
#    par.default <- par(no.readonly=TRUE)
#    if(is.null(colvar) && !is.null(PP$outliers)){
#     colvar <- !PP$outliers
#    }
#    if(is.null(colvar2)){
#     colvar2 <- NULL
#    }
#    par(mfrow=c(2,2), cex =1, cex.lab=1, mai=c(1.2,1.2,1.2,0.4))
#     XLIM0 <- range(PP$CCA$scores[[1]])
#     YLIM0 <- range(PP$CCA$scores[[2]])
#     ni <- min(4, length(PP$CCA$cor))
#     for(i in 1:ni){
#       scatter2D(PP$CCA$scores[[1]][,i], PP$CCA$scores[[2]][,i], 
#                 colvar = colvar, pch = pch,col = col,
#                 xlim = range(c(PP$CCA$scores[[1]][,i], PP$obsCCA[1,i])),
#                 ylim = range(c(PP$CCA$scores[[2]][,i], PP$predCCA[,i])),
#                 xlab="data", ylab="forecast", asp=NA, colkey=FALSE,
#                 main=paste0(i," canonical variables\n(correlation = ", 
#                 round(PP$CCA$cor[i],2),")"))
#       abline(v=PP$obsCCA[1,i], col="chartreuse3", lwd=4)
#       scatter2D(rep(PP$obsCCA[1,i], length( PP$predCCA[,i])), PP$predCCA[,i],
#                 pch=pch2, col=col2, cex=1.5, add = TRUE, colvar = colvar2,
#                 colkey = TRUE)
#       legend("topleft", c("prior","posterior","obs. data"), 
#             col=c("dodgerblue2","firebrick3","chartreuse3"), 
#           pch=c(pch, pch2, NA), lwd = c(NA,NA,2), bty="n", cex=1)
#             grid()
#     }
#     par(par.default)
# }

plotPred <- function(PP, H, Hobs, qlow = 0.05, qhigh = 0.95, thresh = 0.01,
                    ylim = NULL,...,lpos = "topleft"){
par.default <- par(no.readonly=TRUE)

fq005 <- apply(PP$pred,2, function(x) quantile(x, probs = qlow))
fq095 <- apply(PP$pred,2, function(x) quantile(x, probs = qhigh))
fq050 <- apply(PP$pred,2, median)
q005 <- apply(H,2, function(x) quantile(x, probs = qlow))
q095 <- apply(H,2, function(x) quantile(x, probs = qhigh))
q050 <- apply(H,2, median)

# CairoPNG(filename = file.path(paste0(dirProj, "_00_", what, 
#           "_bacteriaPrediction_threshold.png")), 
#          width = 900, height = 900, pointsize = 18*2)
  par(mfrow=c(1,1), mai=c(2.6,2.6,0.4,0.1))
#   if(is.null(ylim)){
    ylim <- range(rbind(H, PP$pred))
#   }
   matplot(t(H),type="l", lty=1, xlab="time (days)", ylim = ylim,
              ylab="Bacteria concentration (g/L)", col="grey", ...)
  lines(fq050, col="firebrick3", lty=1, ...)
  lines(fq095, col="firebrick3", lty="9919", ...)
  lines(fq005, col="firebrick3", lty="9919", ...)
  lines(q050, col="dodgerblue2", lty=1, ...)
  lines(q095, col="dodgerblue2", lty="9919", ...)
  lines(q005, col="dodgerblue2", lty="9919", ...)
  lines(Hobs, col="chartreuse3", ...)
  if(!is.null(thresh)){
    abline(h=thresh, col="black", lty="18", lwd=6)
  }
  legend("topleft", c("prior","posterior", "threshold", "obs. data"), 
         col=c("dodgerblue2","firebrick3", "black", "chartreuse3"), 
         lwd=c(6,6,6,6), 
         lty=c("11","11","18","9919"),bty="n")
	par(par.default)
}

plotDecision <- function(PP, H, Hobs, thresh=0.01){
par.default <- par(no.readonly=TRUE)

exceedCW <- colSums(H > thresh)/nrow(H)
exceedforCw <- colSums(PP$pred > thresh)/nrow(PP$pred)
exceedObs <- which(Hobs > thresh)

# CairoPNG(filename = file.path(paste0(dirProj,"_00_", what, 
#           "_bacteriaContamination.png")), 
#          width = 900, height = 900, pointsize = 18*2)
  par(mfrow=c(1,1), mai=c(1.2,2.8,0.2,0.2))
  plot(exceedCW, type="n", col="dodgerblue2", lwd=2, ylim=c(0,1),
      xlab="time (days)",ylab="contamination probability")
  grid()
  lines(  exceedCW, col="dodgerblue2", lwd=10)
  lines(  exceedforCw, col="firebrick3", lwd=10)
  points(exceedObs, rep(0, length(exceedObs)), pch="*", cex=4, 
          col="chartreuse3")
 legend("topleft", c("prior","posterior","obs. contamination"), 
        col=c("dodgerblue2","firebrick3","chartreuse3"), 
        lwd=c(10,10,NA), pch=c(NA,NA,"*"), pt.cex=c(1,1,4), bty="n")
# dev.off() 
par(par.default)
}

multiFPCA <- function(x, nH=15, nbasis = 70, norder = 3, centerfns = TRUE,
                      scaling = c("none", "global", "local")){
  # 1. basis functions (splines)
  v <- c(1,ncol(x[[1]]))-0.5
  basis <- create.bspline.basis(rangeval=v, nbasis=nbasis, norder)
  x_fd <- multiFdata(x, basisObj = basis)
  if(scaling == "global"){
    x_fd <- fdScale(x_fd, global = TRUE)
  }else if(scaling == "local"){
    x_fd <- fdScale(x_fd, global = FALSE)
  }
  x_fpca <- lapply(x_fd, pca.fd, nharm = nH, centerfns = centerfns)
  if(scaling != "none"){
    for(i in seq_along(x_fpca)){
      x_fpca[[i]]$scale <- x_fd[[i]]$scale
    }
  }
  return(x_fpca)
}

permCCA <- function(X,Y, nperm=100){
#   CCA <- cc(X, Y)
  rho0 <- cc(X, Y)$cor
  n <- nrow(X)
  idx0 <- seq_len(n)
  pp <- numeric(length(rho0))
  rho <- matrix(nrow=nperm, ncol= length(rho0))
  for( i in seq_len(nperm)){
    idx <- sample(idx0, size = n, replace = FALSE) 
    pp <- pp + ( cc(X, Y[idx,])$cor > rho0)
    rho[i,] <- cc(X, Y[idx,])$cor
  }
  # we can improve the estimate of a population proportion by adding 
  # two successes and two failures to the sample.
  pvalue <- (pp + 1)/(nperm + 1)
  return(pvalue)
}

#p.perm(X,Y, rhostart=2)


#' Convert a list of data into a list of functional data
#'
#' @para basisObj a fd basis object
#' @para nH number of harmonic for the fPCA
#' @para x list of matrix (each element of the list = data at one well)
multiFdata <- function(x, basisObj){
  fData <- vector(mode="list", length=length(x))
  vt <- 1:ncol(x[[1]]) - 0.5
  for(i in seq_along(x)){
    fData[[i]] <- smooth.basis(vt, t(x[[i]]), basisObj)$fd
   # DnewSim <- eval.fd(vt, fData[[i]])
   # DnewRes <- sum(as.vector((DnewSim - t(x[[i]]))^2))/prod(dim(DnewSim))
   # cat(i, " : ",DnewRes,"\n")
#     fPCA[[i]] <- pca.fd(fdobj = smbs, nharm = nH, centerfns = TRUE)
  }
  return(fData)
}

fdScale <- function(x, global = TRUE){
  sdfd <- lapply(x, sd.fd)
  glsdfd <- sapply(sdfd, function(x){ mean(x$coefs)})
  if(global){
    glsdfd <- rep(sqrt(sum(glsdfd^2)), length(glsdfd))
  }
  for(i in seq_along(x)){
    x[[i]]$coefs <-  x[[i]]$coefs/glsdfd[i]
    x[[i]]$scale <- glsdfd[i]
  }
  return(x)
}



#' Compute fPCA on a list of matrix
#'
#' @para basisObj a fd basis object
#' @para nH number of harmonic for the fPCA
#' @para x list of matrix (each element of the list = data at one well)
multifPCA <- function(x, nH=15, basisObj){
  fPCA <- vector(mode="list", length=length(x))
  vt <- 1:ncol(x[[1]]) - 0.5
  for(i in seq_along(x)){
    smbs <- smooth.basis(vt, t(x[[i]]), basisObj)$fd
    DnewSim <- eval.fd(vt, smbs)
    DnewRes <- sum(as.vector((DnewSim - t(x[[i]]))^2))/prod(dim(DnewSim))
    cat(i, " : ",DnewRes,"\n")
    fPCA[[i]] <- pca.fd(fdobj = smbs, nharm = nH, centerfns = TRUE)
  }
  return(fPCA)
}


#' Project a vector onto the fPCA space
projfPCA <- function(x, fPCA, nH = NULL){
  xfPCA <- list()
  v <- fPCA[[1]]$meanfd$basis$rangeval
  nbasis <- fPCA[[1]]$meanfd$basis$nbasis
  norder0 <- unlist(strsplit(fPCA[[1]]$meanfd$basis$names[1], '[.]'))[1]
  norder <- as.numeric(substr(norder0, nchar(norder0), nchar(norder0)+1))
  basis <- create.bspline.basis(rangeval=v, nbasis=nbasis, norder)
  for(i in seq_along(x)){
    if(is.null(dim(x[[i]]))){
      vt <- seq_along(x[[i]]) - 0.5
      xx <- (x[[i]])
    }else{
      vt <- 1:ncol(x[[i]]) - 0.5
      xx <- t(x[[i]])
    }
    smbs <- smooth.basis(vt, xx, basis)$fd
    bb <- fPCA[[i]]$harmonics$coefs
    if(!is.null(fPCA[[i]]$scale)){
      smbs$coefs <- smbs$coefs/fPCA[[i]]$scale
    }
    fpcaScores <-  inprod(center.fd2(smbs, fPCA[[i]]$meanfd), 
                          fPCA[[i]]$harmonics)
    if(is.null(nH)){
      nH <- ncol(fpcaScores)
    }
    xfPCA[[i]] <- fpcaScores[,1:nH]
  }
  return(xfPCA)
}

#--- reconstruction artisanale fPCA
# i <- 1
# vt <- 1:ncol(A) - 0.5
# smbs <- smooth.basis(vt, t(A), basisObj)$fd
# ck <- smbs$coef
# PHIbasis0 <-   eval.basis(vt, basisObj)
# x_approx <- PHIbasis0  %*% ck
# sum(abs(t(x_approx) - A))/prod(dim(x_approx))


# # 
# # #' Projection of a vector on PCA-space and backprojection
# # #' @param x numeric vector of length equal to the number of PCA component
# # #' @param PCA list (output from prcomp)
# # #' @param inverse length-one boolean vector. If TRUE the backprojection is 
# # #' computed
# # #' @param comp a integer vector indicating the components selectioned for 
# the
# # #' backprojection
# # projPCA <- function(x, PCA, inverse = FALSE, comp = NULL){
# # #   xNew <- scale(x, PCA$center, PCA$scale)
# # #   t(xNew) %*% PCA$rotation
# #  if(is.null(dim(x))){
# #     dim(x) <- c(1,length(x))
# #   }
# #   if(isTRUE(inverse)){
# #     ROT <- PCA$rotation
# #     if(is.null(comp)){
# #       comp <- 1:min(ncol(ROT),ncol(x))
# #     }
# #     y0 <- x %*%  t(ROT[, comp])
# #     if(!identical(PCA$scale, FALSE)){
# #       y0 <- scale(y0, center = FALSE , scale=1/PCA$scale)     
# #     }
# #     if(!identical(PCA$center , FALSE)){
# #       y0 <- scale(y0, center = -1 * PCA$center, scale=FALSE)    
# #     }
# #     return(y0)
# #   }else{
# #     if(!identical(PCA$center , FALSE)){
# #       x <- scale(x, center = PCA$center, scale=FALSE)    
# #     }
# #     if(!identical(PCA$scale, FALSE)){
# #       x <- scale(x, center = FALSE , scale=PCA$scale)     
# #     }
# #     x %*%  PCA$rotation
# # #     xCentered <- (x - PCA$center)
# # #     if(identical(fDPCA$scale, FALSE)){
# # #       t( xCentered) %*% PCA$rotation
# # #     }else{
# # #       t(xCentered * fDPCA$scale) %*% PCA$rotation
# # #     }
# #   }
# # }
# # 
# # projCCA <- function(x, CCA, type=c("x","y"), inverse = FALSE, 
# #                     mu = 0.00001, center = NULL){
# # 
# #   type <- match.arg(type, c("x","y"))
# #   if(isTRUE(inverse)){
# #     if(type == "y"){
# #       if(is.null(dim(x))){
# #         dim(x) <- c(length(x),1)
# #       }
# #       B <- CCA$ycoef
# #       if(is.null(dim(B))){
# #         dim(B) <- c(length(B),1)
# #       }
# #       B <- B[,1:ncol(x), drop=FALSE]
# #       if(!is.null(center)){
# #         x <- scale(x, center = -1*center[1:ncol(x)], scale = FALSE)
# #       }else{
# #       }
# #       y <- try( x %*% t(B) %*% solve(B %*% t(B)), silent =TRUE)
# #       if(class(y) == "try-error"){
# #         y <- x %*% t(B) %*% 
# #               solve(B %*% t(B) + diag(mu, nrow(B)))
# #       } 
# #       return(y)
# #     }else{
# #       stop("not yet implemented\n")
# #     }
# #   }else{
# #     if(type == "x"){
# #       if(is.null(dim(x))){
# #         dim(x) <- c(1, length(x))
# #       }
# #       if(!is.null(center)){
# #         x <- scale(x, center = center, scale = FALSE)
# #       }
# #       x[, 1:nrow(CCA$xcoef)] %*% CCA$xcoef
# #     }else{
# #       stop("not yet implemented\n")
# #     }
# #   }
# # }

fPCAScores <- function(fPCA, nH = NULL){
  fD <- c()
  if(is.null(nH)){
    nH <- ncol(fPCA[[1]]$scores)
  }
  for(i in seq_along(fPCA)){
    fD <- cbind(fD, fPCA[[i]]$scores[,1:nH])
  }
  return(fD)
}


levFilter <- function(X, fac = 2){
  Dc_l <- apply(Dc, 2, leverage)
    test <- apply(Dc_l < fac*2/nrow(Dc),1,any)
    return(test)
}

leverage <- function(x){
  xmean <- mean(x)
  1/length(x) + ((x - xmean)^2) / (sum((x-xmean)^2))
}
stdRes <- function(x,y, intrcp = FALSE){
  if(isTRUE(intrcp)){
    xx <- cbind(rep(1,length(x)), x)
  }else{
    dim(x) <- c(length(x),1)
    xx <- x
  }
  dim(y) <- c(length(y),1)
  #G <- solve(t(y) %*% y + 0.0000001 , t(y) %*% x)
  G <- solve(t(xx) %*% xx  , t(xx) %*% y)
  res <- y - xx %*% G
  s <- sqrt(1/(length(x)-2) * sum(res^2))
  return(res/(s * sqrt(1 - xl)))
}

lmf <- function(x, y, levTsh=2, stdResTsh = 3.5){
  # leverage points
  xl <- leverage(x)
  test <- xl > levTsh*2/length(x)
  # standardise residuals
  # solve Y = G X for G (ax + b = y)
  dim(y) <- c(length(y),1)
  dim(x) <- c(length(x),1)
  xx <- cbind(rep(1,length(x)), x)
  xx <- x
  #G <- solve(t(y) %*% y + 0.0000001 , t(y) %*% x)
  G <- solve(t(xx) %*% xx  , t(xx) %*% y)
  res <- y - xx %*% G
  s <- sqrt(1/(length(x)-2) * sum(res^2))
  res2 <- res/(s * sqrt(1 - xl))
  test2 <- abs(res2) > stdResTsh
  test3 <- abs(res2) > 4* sd(res2)
  test_all <- !((test2 & test) | test3)
  return(test_all)
}

lmRes <- function(x,y,intrcp = FALSE){
if(isTRUE(intrcp)){
    xx <- cbind(rep(1,length(x)), x)
  }else{
    dim(x) <- c(length(x),1)
    xx <- x
  }
  dim(y) <- c(length(y),1)
  #G <- solve(t(y) %*% y + 0.0000001 , t(y) %*% x)
  G <- solve(t(xx) %*% xx  , t(xx) %*% y)
  res <- y - xx %*% G
  return(res)
}


lmfwrap <- function(i, x, y, levTsh=2, stdResTsh = 3.5){
  return(lmf(x[,i], y[,i], levTsh = levTsh, stdResTsh= stdResTsh))

}


# # CHAP 14 from Bayesian Data Analysis
# # Third Edition
# # ISBN 978-1-4398-9820-8
# # bayesian linear regression y = ax + b
# blrfit <- function(y, x ,N = 10){
#   lmfit <- lm(y~x)
#   QR<-lmfit$qr
#   df.residual<-lmfit$df.residual
#   R<-qr.R(QR) ## R component
#   coef<-lmfit$coef
#   Vb<-chol2inv(R) ## variance(unscaled)
#   s2<-(t(lmfit$residuals)%*%lmfit$residuals)
#   s2<-s2[1,1]/df.residual
# 
#   ## now to sample residual variance
#   sigma <- df.residual*s2/rchisq(N, df.residual)
#   coef.sim <-sapply(sigma, function(x)  MASS::mvrnorm(1, coef, Vb*x))
#   ret <- data.frame(t(coef.sim))
#   names(ret) <- names(lmfit$coef)
#   ret$sigma <- sqrt(sigma)
#   return(ret)
# }
# blrpred <- function(y, x , xobs, N = 1000){
#   lmfit <- lm(y~x)
#   QR <- lmfit$qr
#   df.residual <- lmfit$df.residual
#   R <- qr.R(QR) ## R component
#   coef <- lmfit$coef
#   Vb <- chol2inv(R) ## variance(unscaled)
#   s2<-(t(lmfit$residuals)%*%lmfit$residuals)
#   s2<-s2[1,1]/df.residual
# 
#   ## now to sample residual variance
#   sigma2 <- df.residual*s2/rchisq(N, df.residual)
#   coefsim <-sapply(sigma2, function(x)  MASS::mvrnorm(1, coef, Vb*x))
#   # ret <- data.frame(t(coef.sim))
#   # names(ret) <- names(lmfit$coef)
#   # ret$sigma <- sqrt(sigma2)
#   #if(is.null(xobs)){
#   #  dim(xobs) <- c(length(xobs),1)
#   #}
#   xo <- matrix(1, nrow = 1, ncol=length(xobs) + 1)
#   xo[,2:ncol(xo)] <- xobs
#   dd <- xo %*% coef
#   ss <- s2 * (1 + xo %*% Vb %*% t(xo))
#   dof <- length(x) - length(coef)
#   ypred <- mvtnorm::rmvt(n = N, delta = dd, sigma = ss, df = dof)
#   return(ypred)
# #   xo <- matrix(1, nrow = N, ncol=length(xobs) + 1)
# #   xo[,2:ncol(xo)] <- xobs
# #   mu <- xo %*% coefsim
# #   V <- (1 + xobs %*% Vb %*% t(xobs)) * ret[,ncol(ret), drop = FALSE]
# }


#### OLD
# posterior of multi-variate regression with Gaussian prior
# and Gaussian modeling error
# BUT HERE, THE ERROR ARE UNCORRELATED, THE Y AND X ARE UNCORRELATED
# if xobs is multi-dimensional: m x n, with m = number of observations
#   return a multi-dimensional mean of the same size (m x n)
postGaussian <- function(X, Y, xobs, outliers=TRUE, levTsh = 2, 
                        stdResTsh = 3.5){
  if(is.null(dim(xobs))){
    dim(xobs) <- c(1, length(xobs))
  }
  if(isTRUE(outliers)){
    X_l <- sapply(1:ncol(X), lmfwrap, X, Y, levTsh = levTsh, 
                  stdResTsh = stdResTsh)
    test <- apply(X_l,1,all)
    X <- X[test,, drop = FALSE]
    Y <- Y[test,, drop = FALSE]
  }
  C_Y <- cov(Y)
  if(ncol(Y) == 1){
    C_Y <- diag(ncol(X))
  }
  Y_Mean <- matrix(colMeans(Y), nrow=ncol(Y),ncol=nrow(xobs))
#   dim(Y_Mean) <- c(length(Y_Mean),nrow(xobs))

  # Find best linear bit between X and Y
  # D = G H
  G <- solve(t(Y) %*% Y + diag(0.0000001, ncol(Y)) , t(Y) %*% X)
  XSim <- Y %*% G
  # sum(abs((XSim) - X))

  XDiff <- (X - (XSim))
  C_T <- t(XDiff) %*% XDiff/nrow(X)

  mu_post <- Y_Mean + C_Y %*% t(G) %*% 
                pinv(G %*% C_Y %*% t(G) + C_T) %*% 
                (t(xobs) - G %*% Y_Mean)

  C_post <- inv(t(G) %*% pinv(C_T) %*% G + inv(C_Y))
  
  return(list("mean" = t(mu_post), "cov" = C_post, "outliers"=test))
}


# posterior of multi-variate regression with Gaussian prior
# and Gaussian modeling error
# BUT HERE, THE ERROR ARE UNCORRELATED, THE Y AND X ARE UNCORRELATED
# if xobs is multi-dimensional: m x n, with m = number of observations
#   return a multi-dimensional mean of the same size (m x n)
lmpred <- function(X, Y, xobs, npred = 10, outliers=TRUE, levTsh = 2, 
                        stdResTsh = 3.5, w = 0.25, localw = FALSE){
  if(is.null(dim(xobs))){
    dim(xobs) <- c(1, length(xobs))
  }
  if(isTRUE(outliers)){
    X_l <- sapply(1:ncol(X), lmfwrap, X, Y, levTsh = levTsh, 
                  stdResTsh = stdResTsh)
    test <- apply(X_l,1,all)
    X <- X[test,, drop = FALSE]
    Y <- Y[test,, drop = FALSE]
  }
  if(localw == TRUE){
    i <- 1
    yi <- Y[,i]
    yimean <- mean(yi)
    xi <- X[,i]
    lmfit <- lm(yi ~ xi)
    QR <- lmfit$qr
    df.residual <- lmfit$df.residual
    R <- qr.R(QR) ## R component
    coef <- lmfit$coef
    Vb <- chol2inv(R) ## variance(unscaled)
    s2 <- (t(lmfit$residuals) %*% lmfit$residuals)
    s2 <- s2[1,1]/df.residual
    bf<-bayesfit(lmfit, 10)
  }
  C_Y <- cov(Y)
  Y_Mean <- matrix(colMeans(Y), nrow=ncol(Y),ncol=nrow(xobs))
  # Find best linear bit between X and Y
  # D = G H
  G <- solve(t(Y) %*% Y + diag(0.0000001, ncol(Y)) , t(Y) %*% X)
  XSim <- Y %*% G
  # sum(abs((XSim) - X))

  XDiff <- (X - (XSim))
  C_T = t(XDiff) %*% XDiff/nrow(X)

  mu_post <- Y_Mean + C_Y %*% t(G) %*% 
                pinv(G %*% C_Y %*% t(G) + C_T) %*% 
                (t(xobs) - G %*% Y_Mean)

  C_post <- inv(t(G) %*% pinv(C_T) %*% G + inv(C_Y))
  Xpred <- MASS::mvrnorm(n = npred, mu = t(mu_post), Sigma = C_post)
  return(list("pred" = Xpred,  "outliers"=test))
}




# center a functional object "fobj" by another functional object "fdobj2"
# (useful when "fdobj2" represents the mean of "fdobj" and you want to 
# remove this mean from "fdobj").
center.fd2 <- function (fdobj, fdobj2) 
{
    coef <- as.array(fdobj$coefs)
    coef2 <- as.array(fdobj2$coefs)
    coefd <- dim(coef)
    ndim <- length(coefd)
    basisobj <- fdobj$basis
    nbasis <- basisobj$nbasis
    if (ndim == 2) {
        coefmean <- coef2
        coef <- sweep(coef, 1, coefmean)
    }
    else {
        nvar <- coefd[3]
        for (j in 1:nvar) {
            coefmean <- apply(coef[, , j], 1, mean)
            coef[, , j] <- sweep(coef[, , j], 1, coefmean)
        }
    }
    fdnames <- fdobj$fdnames
    fdnames[[3]] <- paste("Centered", fdnames[[3]])
    centerfdobj <- fd(coef, basisobj, fdnames)
    return(centerfdobj)
}

#--- normal score transformation
# cf. technical note: the normal quantile transformation and its application
# in a flood forecasting system.
# check tRank {multic}
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2367615/
qTrans <- function(x){
  # the line below is equivalent to
  # 1. sort the samples from smallest to largest > x <- sort (x)
  # 2. estimate cumulative probability using a plotting position like
  #     i /(n+1)  >  y <- ppoints(x, a=0)
  y2 <- qqnorm(x, plot.it = FALSE)$x
  y_n <- approx(x, xout = y2, rule = 2)$y
  
}


# histogram transformation to normal distributed data with same mean and sd
# https://msu.edu/~ashton/temp/nscore.R
nScoreTrans <- function(x, inverse = FALSE, tbl = NULL){
  if(isTRUE(inverse)){
    if(is.null(tbl)){
      stop("tbl must be provided")
    }
    min_x <- min(tbl[,1])
    max_x <- max(tbl[,1])
    min_sc <- min(x)
    max_sc <- max(x)
    x1 <- c(min_x, tbl[,1], max_x)
    nsc <- c(min_sc, tbl[,2], max_sc)
    
    back.xf <- approxfun(nsc,x1) # Develop the back transform function
    val <- back.xf(x)
    return(val)
  }else{
#     fx <- ecdf(x)
#     x1 <- head(knots(fx), -1)
#     xCDF <- fx(x1)
#     y <- qnorm(xCDF, mean(x), sd = sd(x))
    if(!is.null(tbl)){
      x1 <- tbl[,1]
      y  <- tbl[,2]
    }else{
      x1 <- sort(x)
      y <- sort(qqnorm(x, plot.it = FALSE)$x)
      tbl <- data.frame(x = x1, nscore = y)
    }
    y_n <- approx(x1, y, xout = x, rule = 2)$y
    return(list(x = y_n, table = tbl))
  }
}



dodGSAForMe <- function(D, P, ncl = 5, dirProj, FUN, pfx="DATA_"){
  
  #--- distance
  dData <- dist(D)

  #---- clustering (PAM)
  kmed3 <- cluster::pam(dData, ncl)
  pc <- kmed3$clustering

  #--- MDS
  dataMDS <- cmdscale(dData, k = 13, eig = TRUE)
  CairoPNG(filename = file.path(dirProj, 
          paste0(pfx,"MDS_cl_",ncl,".png")), 
          width = 1280, height = 660, pointsize = 14)
      par(mfrow=c(1,3))
      plot((cumsum(dataMDS$eig)/sum(dataMDS$eig))[1:20], type="o", pch=20,
            ylim=c(0,1), xlab="eigenvalues", ylab="cumulative variance")
      points2D(dataMDS$points[,1],dataMDS$points[,2],pch=20,colvar=pc,asp=1)
      mtext(paste0("Coloured by cluster group"),side=3,at=par("usr")[1],
            line=1, adj=c(0,0), cex=1.2, col="black")
      points3D(dataMDS$points[,1],dataMDS$points[,2],dataMDS$points[,3],
              colvar=pc,      bty="f",cex=1,pch=20, clab = "cluster",
              ticktype="detailed",       theta=40 , expand=1,scale=FALSE,
              xlab="axe 1",ylab="axe 2", zlab="axe 3",shade=TRUE,
              border="black", colkey = list(width = 0.5, length = 0.5,
              cex.axis = 0.8, side = 1), col.axis="black",col.panel="white",
              col.grid="grey",lwd.panel=1,lwd.grid=2,box=TRUE)
  dev.off()

  #--- dGSA 
  # main factor + p-values
  MF <- apply(P,2,distCDF, cl=pc, nPerm = 1000, n=99)
  
  CairoPDF(file = file.path(dirProj, paste0(pfx,"dGSA_",ncl,".png")),
          width = 8, height = 16, pointsize = 14)
    dGSABarplot(MF)
  dev.off()
  CairoPNG(filename = file.path(dirProj, paste0(pfx,"dGSA_clusts_",ncl,".png")),
          width = 680, height = 960, pointsize = 14)
    GSABarplotClust(MF)
  dev.off()
  
#   mf <- distCDF(P[,1], , cl=pc, nPerm = 100, n=99)
#   
#   #---- permutation test
#   cluster <- parallel::makeCluster(4, type = "SOCK")
#   doSNOW::registerDoSNOW(cluster)
#   cl <- pc
#   n <- 99
#   nPerm <- 1000
#   alpha <- c(0.9,0.95,1)
#   x0 <- MF[,i]
#   snow::clusterExport(cluster, c("cl", "n", "nPerm","alpha", "x0"), 
#                       envir = environment())
#   permMF <- foreach(i=1:ncol(P), .combine=cbind) %dopar% {
#     FUN(P[,i], MF[,i], cl=cl,n=n,nPerm=nPerm)
#   }
#   stopCluster(cluster)
# 
#   # normalisation (sensitivity main factors)
#   MF090 <- MF/permMF[,seq(1,length.out=ncol(P),by=3)]
#   MF095 <- MF/permMF[,seq(2,length.out=ncol(P),by=3)]
#   MF100 <- MF/permMF[,seq(3,length.out=ncol(P),by=3)]


#   CairoPNG(filename = file.path(dirProj, paste0(pfx,"dGSA_",ncl,".png")),
#           width = 680, height = 960, pointsize = 14)
#     GSABarplotClust(MF095, MF090, MF100, 
#                     main=paste0("data (",ncl, " clusters)"))
#   dev.off()

}


ifft <- function(x){
  return( Re(fft(x, inverse=TRUE))/length(x) )
}


