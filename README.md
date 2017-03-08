# gwModBac

Groundwater flow simulation and particle tracking to forecast microbial concentration in a drinking water extraction well.

## Installation

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


