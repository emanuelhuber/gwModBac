# gwModBac

Groundwater flow simulation and particle tracking to forecast microbial concentration in a drinking water extraction well.

## Installation/Run

### Requirement
- need to have "gdal_polygonize.py" installed (--> `install python-gdal`)
- need to have MODFLOW and MODPATH installed (`mfusg` for MODFLOW usgs,  `mp6` for MODPATH). If you need binaries compiled for linux (ubuntu) please contact me
- need to have [https://cran.r-project.org/](R) installed

### Run
Several possiblities:

* Open R, copy-paste `run.R` or source `run.R` (but don't forget to adapt line 15 in `run.R` the `DIR` variable to your directory structure
* Open terminal, enter 
    `Rscript --vanilla ~/simulations/gwModBac/run.R test.txt 50 15 2 1`
    * `~/simulations/gwModBac/run.R` --> filepath of `run.R`
    * `test.txt` path of the output file (here, `test.txt` will be create in the directory `~/simulations/gwModBac/`
    * `50` --> nx: number of cells along the x-direction
    * `15` --> ny: number of cells along the y-direction
    * `2` --> nz: number of cells along the z-direction
    * `1` --> number of runs (but one output file with all the output together)
    
### Notes
Recommanded minimal grid size:

* nx >= 45
* ny >= 15
* nz >= 2

## What does "run.R" Do?
1. Read the hyperparameters/parameters from `para.R`
2. Read the data (river stage time-series, groundwater head time-series) from `data/`
3. Create the simulation grid (according to parameters)
4. Define 
   * the river properties (location, stage, riverbed), 
   * the location of the specified head boundary, 
   * the particles to estimate the microbila concentration in the drinking water extraction well.
5. Run Monte Carlo simulations (unconditional)
	
    i. simulate the hydraulic properties: hydraulic conductivity (Gaussian random field), porosity, etc.
    ii. simulate the specified head boundary conditions (Gaussian process)
        for the 10-days ahead forecast:
		
        a. simulate the precipitation for the next 10 days
        b. using a convolution model between precipitation and river stage, simulate the river stage for the next 10 days
        c. using a convolution model between river stage and groundwater
			heads at the observation wells, simulate the groundwater heads 
			at the observation wells for the next 10 days
        d. Simulate the head boundary conditions conditional on the
			groundwater heads the observation wells
        e. Interpolate from the head boundary conditions the initial
			heads
    iii. run MODFLOW
    iv. run MODPATH to simulate pathway of the microbes that reach the well drinking water extraction well
    v. Simulate the microbial concentration in the river from the river stage
    vi. Predict the microbial concentration in the well according to an exponential decay model that simulates the mechanical filtration of microbes in the aquifer


