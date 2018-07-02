


uni <- c("days","meters")
wetting <- c("wetfct" = 0.1 , "iwetit" = 2  , "ihdwet" = 0, "wetdry"= 0.1)
cst_mps2mpd <- 3600*24
timeFormat <- "%Y%m%d"
nstp_pred <- 10              # day-ahead prediction
nstp <- 21 + nstp_pred       # number of time steps

timeSel <- seq(139, by=1, length=nstp)


# model size + default resolution
modGrid <- list(x  = c(min = 0, max = 50),   # 0 m to 200 m along x axis
                y  = c(min = 0, max = 150),
                z  = 10,
                nx = 25,      # number of cols (x axis)
                ny = 75,      # number of rows (y axis)
                nz = 5)       # number of cells (z axis)

# river
river <- list(perimeter = cbind(c(modGrid$x["min"], 4,   4,   modGrid$x["min"]),
                                c(modGrid$y["min"], modGrid$y["min"], 
                                  modGrid$y["max"], modGrid$y["max"])),
              slope     = 2/100,
              stepL     = 40,
              stepH     = 0.5,
              minStage  = 0.1,
              depth     = 0.5,
              bedT      = 0.25)  # river bed thickness
river$stepS <- (river$stepH - river$slope*river$stepL)/river$stepL
nrivObs <- 10
# number of constraint for derivative of GP
nbc <- 10


#--- Extraction Well (forecast = bacteria concentration)
wellExt <- list(  x = c(20.5) + 5,
                  y = c(25.5)  + 25,
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

#--- convolution parameters (unit hydrograph)
convMod <- list(mu = 0.05,
                nf = 20)

#---- PRIOR > parameters for the probability density functions
prior <- list(h_sig   = list(type = "unif", min = 0.05, max = 0.1),
              h_lx    = list(type = "unif", min = 20, max = 60),
              h_vx    = list(type = "unif", min = 1, max = 3),
              h_hx    = list(type = "unif", min = 0.05, max = 1),
              h_lt    = list(type = "unif", min = 0.25, max = 1.5),
              h_ht    = list(type = "unif",  min = 0.01, max = 0.2),
              #h_scale = list(type = "unif", min = 2, max = 10),
              Cr      = list(type = "norm", mean = -2, sd = 0.75),
              ss      = list(type = "unif", min = 1*10^-6, max = 10*10^-5),
              sy      = list(type = "unif", min = 0.2, max = 0.35),
              p       = list(type = "unif", min = 0.25, max = 0.35),
              K_mean  = list(type = "norm", mean = -3, sd = 0.1),
              K_sd    = list(type = "unif", min = 0.1, max = 0.75),
              K_l     = list(type = "unif", min = 5, max = 15),
              K_nu    = list(type = "unif", min = 1.5, max = 3),
              K_hani  = list(type = "unif", min = 100, max = 150),
              #K_vani  = list(type = "unif", min = 3, max = 7)
              K_hstr  = list(type = "unif", min = 1/20, max = 1/2),
              K_vstr  = list(type = "unif", min = 15, max = 30),
              K_nug   = list(type = "unif", min = 0.05, max = 0.2)
             )

