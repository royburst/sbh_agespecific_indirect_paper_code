## #######################################################
## AUTHOR: Roy Burstein
## DATE: JANUARY 2018
## PURPOSE: makes dated config and launches EXTVAL_run_models_fitting.R
## #######################################################


## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ##
#  MAKE SURE YOU HAVE UPDATED RUN DATE FOR THE BEST MODEL RUN (IF NEEDED)
#  FILE LOCATED HERE: /share/geospatial/mbg/u5m/died/input/sbh/dated_sbh_model_fitting/best_run_date.txt
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ##


# Assumes you have your repo cloned into /share/geospatial/USERNAME/code/
user      <- Sys.info()['user']
gitbranch <- 'master'


# git pull, if you want
system(sprintf('cd /share/code/geospatial/%s/u5m\ngit pull origin %s', user, gitbranch))


# root dir for outputs
root <- '/share/geospatial/mbg/u5m/died/input/sbh/dated_sbh_model_fitting'

# load custom functions
setwd(sprintf('/share/code/geospatial/%s/u5m/bh_prep', user))
source('sbh_utils.R')

# get a run date
run_date <- make_time_stamp()

# make directory for the rundate
make_output_dirs(root = root, run_date = run_date)
rd_dir <- sprintf('%s/%s',root,run_date)



#######################################################
## ADD TO NOTES FILE FOR THIS RUN_DATE
write_rd_notes(note = 'Running predictions with censored fix in place plus svy weights in the model fit. ')



#######################################################
## SET OPTIONS AND WRITE CONFIGS
# CONFIG Formula
zzz_formula <- cbind('formula',
                     'died ~ -1 + s(yrborn, sdi, by = factor(ab), k = 9) + s(cdceb100, birthorder, mothage_atbirth, k = 9) + factor(ab) + s(nid, bs = \'re\') + s(country_ab, bs = \'re\')')

# CONFIG Cores
zzz_cores <- cbind('cores', '10')


# CONFIG base year, typically 1970
zzz_by <- cbind('base_year','1970')


# place options into config and write
config <- matrix(NA,ncol=2,nrow=0)
for(o in ls()[grepl('zzz',ls())]) config <- rbind(config, get(o))
write.csv(config,file=sprintf('%s/config.csv',rd_dir))

## Age bins, seperate config - should match whatever is being prepped for CBH
ab_times <- data.frame(
  tstart = c(0, 1,  6, 12, 24, 36, 48),
  tstop  = c(1, 6, 12, 24, 36, 48, 60),
  ab     = c('NN','PNN1','PNN2','1yr','2yr','3yr','4yr') )

ab_times$ab <- paste0(1:length(ab_times$ab),"_",ab_times$ab)
write.csv(ab_times,file=sprintf('%s/ab_times.csv',rd_dir))



#######################################################
## Submit the script for each regional run

# Loop through regions
regions <- c('NAME','SSAWC','ASIA','SSASE','LAC')

for(region in regions){
  # make qsub
  qsub <- make_qsub(rtdir       = root,
                    rd          = run_date,
                    cores       = as.numeric(zzz_cores[1,2]),
                    memory      = 150,
                    proj        = 'proj_geospatial',
                    singularity = TRUE,
                    coderepo    = sprintf('/share/code/geospatial/%s/u5m', user),
                    codescript  = 'EXTVAL_run_models_fitting.R',
                    geo_nodes   = TRUE,
                    cntry       = region  )

  system(qsub)
}


system('qstat')
