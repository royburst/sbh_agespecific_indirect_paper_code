## #######################################################
## AUTHOR: Roy Burstein
## DATE:  2018
## PURPOSE: launches EXTVAL_run_models_predict.R per NID
## #######################################################

####################################################################
# USER INPUT #######################################################
# IF YOU WANT TO ADD SURVEYS TO AN EXISTING RUN DATE SET IT HERE, OTHERWISE SET IT TO NULL
existing_run_date <- NULL

####################################################################


# Assumes you have your repo cloned into /share/geospatial/USERNAME/code/
user      <- Sys.info()['user']
gitbranch <- 'master'

# git pull, if you want
system(sprintf('cd /share/code/geospatial/%s/u5m\ngit pull origin %s', user, gitbranch))

# root dir for outputs
root <- '/share/geospatial/mbg/u5m/died/input/sbh/dated_sbh_model_prediction'

# load custom functions
setwd(sprintf('/share/code/geospatial/%s/u5m/bh_prep', user))
source('sbh_utils.R')
library(data.table)
library(splitstackshape)

# get a run date
if(is.null(existing_run_date)){
  run_date <- make_time_stamp()
  make_output_dirs(root = root, run_date = run_date)
  message(sprintf('New run date: %s',run_date))
} else {
  run_date <- existing_run_date
  message(sprintf('Using existing run date: %s',run_date))
}

# make directory for the rundate
rd_dir <- sprintf('%s/%s/',root,run_date)

# grab all nids of SBH surveys in our database
nids <- fread('/home/j/WORK/11_geospatial/02_processed data/U5M/Global/data_list.csv')
nids <- subset(nids, type == 'SBH') # keep only those that are SBH
nids <- subset(nids, year >= 1998 ) # keep only those after 1998


#######################################################
# QSUB all NIDS

# Set the below to TRUE  if you want estimates at the survey level (use this for SBH-CBH paper)
# Set the below to FALSE if you want estimates at the smallest available geography level
aggregate_to_nid <- TRUE


for(n in nn){ 
  
  # make qsub
  qsub <- make_qsub(rtdir       = root,
                    rd          = run_date,
                    cores       = 10, 
                    memory      = 150,
                    proj        = 'proj_geospatial',
                    coderepo    = sprintf(
                      '/share/code/geospatial/%s/u5m', user),
                    codescript  = 'EXTVAL_run_models_predict.R',
                    geo_nodes   = TRUE,
                    singularity = TRUE, # set to 10 cores right now
                    cntry       = n,
                    adl_args    = aggregate_to_nid)

  system(qsub)
}

system('qstat')


### NOTE: ONCE THESE MODELS HAVE FINISHED, YOU NEXT RUN EXTVAL_SBH_COLLAPSE.R TO COMBINE THEM
