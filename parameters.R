


uni <- c("days","meters")
wetting <- c("wetfct" = 0.1 , "iwetit" = 2  , "ihdwet" = 0, "wetdry"= 0.1)
cst_mps2mpd <- 3600*24
timeFormat <- "%Y%m%d"
nstp_pred <- 10              # day-ahead prediction
nstp <- 21 + nstp_pred       # number of time steps
timeSel <- seq(139, by=1, length=nstp)
modGrid <- list(L  = c(min = 0, max = 150),   # 0 m to 200 m along x axis
                W  = c(min = 0, max = 50),
                H  = 10,
                nx = 75,      # number of cells (x axis)
                ny = 25,       # number of cells (y axis)
                nz = 5)       # number of cells (z axis)
river <- list(perimeter = cbind(c(modGrid$W["min"], 4,   4,   modGrid$W["min"]),
                                c(modGrid$L["min"], modGrid$L["min"], 
                                  modGrid$L["max"], modGrid$L["max"])),
              slope     = 2/100,
              stepL     = 40,
              stepH     = 0.5,
              minStage  = 0.1,
              depth     = 0.5,
              bedT      = 0.25)  # river bed thickness
river$stepS <- (river$stepH - river$slope*river$stepL)/river$stepL
# piezometers > hydraulic head observations
piez <- list(x = cbind(c(17,   44.8,  9.1, 30.5, 13.4, 22.9, 36.2),
                       c(19.8, 12.5, 35.6, 44.5, 58.6, 72.5, 68) + 15,
                       c(14.9, 15.2, 14.7, 14.4, 14.6, 15.4, 16.0)))
nrivObs <- 10
nbc <- 10
#--- convolution parameters (unit hydrograph)
convMod <- list(mu = 0.05,
                nf = 20)
#--- Extraction Well (forecast = bacteria concentration)
wellExt <- list(  x = c(20.5) + 5,
                  y = c(25.5) + 25,
               ztop = c(16.2),
               zbot = c(14),
                 id = c(1),
                  q = c(-1000))
#---- particles forecast ----#
wellExtPart <- list(nxy = 25,  # number of particles horizontally
                   nz  = 10,  # number of "layers" of particles vertically
                   # particle time release
                   t   = nstp - (nstp-nstp_pred):(nstp -1))  
                   
bactConcPara <- list(lambda = 0.1,
                     a      = 1,
                     b      = 0.06,
                     Cpmin  = 0)

#---- PRIOR > parameters for the probability density functions
prior <- list(h_sig   = list(type = "unif", min = 0.05, max = 0.1),
              h_lx    = list(type = "unif", min = 20, max = 60),
              h_vx    = list(type = "unif", min = 1, max = 3),
              h_hx    = list(type = "unif", min = 0.05, max = 1),
              h_lt    = list(type = "unif", min = 0.25, max = 1.5),
              h_ht    = list(type = "unif",  min = 0.01, max = 0.2),
              h_scale = list(type = "unif", min = 2, max = 10),
              Cr      = list(type = "unif", mean = -2, sd = 0.75),
              ss      = list(type = "unif", min = 1*10^-6, max = 10*10^-5),
              sy      = list(type = "unif", min = 0.2, max = 0.35),
              p       = list(type = "unif", min = 0.25, max = 0.35),
              K_mean  = list(type = "norm", mean = -3, sd = 0.1),
              K_sd    = list(type = "unif", min = 0.1, max = 0.75),
              K_l     = list(type = "unif", min = 5, max = 15),
              K_nu    = list(type = "unif", min = 1.5, max = 3),
              K_hani  = list(type = "unif", min = 100, max = 150),
              K_hstr  = list(type = "unif", min = 1/20, max = 1/2),
              K_vstr  = list(type = "unif", min = 15, max = 30),
              K_nug   = list(type = "unif", min = 0.05, max = 0.2),
              K_vani  = list(type = "unif", min = 3, max = 7))

