# Code for *PLOS Medicine* Article "Development and validation of a new method for indirect estimation of neonatal, infant, and child mortality" 

This repository hold all code used to run the analysis and produce the figures for this paper. 

## Basic information

Scripts in this repo begin with either `CROSSVAL`, `EXTVAL` or `MANUSCRIPT`. Scripts with the prefix `CROSSVAL` were used to run and fit the models which were used in the cross-validation component of the analysis, while the `EXTVAL` scripts were used for the external validation. `MANUSCRIPT` scripts are where the figures get made. 

Note that much of the code here was written to run specifically in the Univa Grid Engine environment for parallel computing. As such several of the larger processes have two scripts, one noted in the names `launch`, and another `run`. The `launch` script gives the UGE scheduler instructions on what data to run the process, typically looping over surveys or countries, while the `run` script is the general processing instructions for each chunk of data. 

Many of the functions used throughout the analysis will be found in `utils.R`


