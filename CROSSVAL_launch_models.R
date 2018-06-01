## Roy Burstein
## This script is used to launch many parallel jobs of CROSSVAL_holdout_run_parallel.R


system('cd /share/geospatial/royburst/code/sbh_cbh\ngit pull origin master') # pulls from develo

# This launches validation models en masse

# load custom functions
setwd('/share/geospatial/royburst/code/sbh_cbh')
source('utils.R')

# root dir for outputs
root <- '/share/geospatial/royburst/sbh_cbh/output'

# get a run date
run_date <- make_time_stamp()

# make directory for the rundate
make_output_dirs(root = root, run_date = run_date)
rd_dir <- sprintf('%s/%s',root,run_date)




#######################################################
#####  ADD TO NOTES FILE FOR THIS RUN_DATE
write_rd_notes(note = 'INDIV with censfix')


#######################################################
##### SET OPTIONS AND WRITE CONFIGS
## CONFIG Formula

# FULL
#zzz_formula <- cbind('formula',
#  'died ~ -1 + s(yrborn, sdi, by = factor(ab), k = 9) + s(cdceb100, birthorder, mothage_atbirth, k = 9) + factor(ab) + s(nid, bs = \'re\') + s(country_ab, bs = \'re\')')

# INT
#zzz_formula <- cbind('formula',
#                     'died ~  -1 + yrborn*sdi*factor(ab) + cdceb100*birthorder*mothage_atbirth + s(nid, bs = \'re\') + s(country_ab, bs = \'re\')')

# ADD
#zzz_formula <- cbind('formula',
#                     'died ~  -1 + yrborn + sdi + factor(ab) + cdceb100 + birthorder + mothage_atbirth + s(nid, bs = \'re\') + s(country_ab, bs = \'re\')')

# TREND
#zzz_formula <- cbind('formula',
#                     'died ~  -1 + s(yrborn, sdi, by = factor(ab), k = 9) + factor(ab) + s(nid, bs = \'re\') + s(country_ab, bs = \'re\')')

# INDIV
zzz_formula <- cbind('formula',
                     'died ~  -1 + s(cdceb100, birthorder, mothage_atbirth, k = 9) + factor(ab) + s(nid, bs = \'re\') + s(country_ab, bs = \'re\')')







## CONFIG Cores
zzz_cores <- cbind('cores', '15')

## CONFIG training extent: 'survey', 'country', 'region', 'global'
  # NOTE if te is country, remove country or country_ab random effect from the formula
zzz_te <- cbind('training_extent', 'region')

## CONFIG predict for, oos or is or both 'c(\'oos\',\'is\')'
  # NOTE zzz_te = 'survey' and predict_for = 'is' gives a one survey run with in sample predictions
  # Typically we are interested in oos fits, likely on a regional training extent
zzz_pf <- cbind('predict_for','c(\'oos\')')

## CONFIG base year, typically 1970
zzz_by <- cbind('base_year','1970')

## CONFIG hold out strategy, only last year is supported
zzz_ho <- cbind('ho_strat','last_year')

## CONFIG number of draws for prediction
zzz_nd <- cbind('ndraws','1000')

## CONFIG number of years with which to aggregate over estimates. 1 is no agg
zzz_ya <- cbind('yrs_to_agg','1')


## CFG END place options into config and write
config <- matrix(NA,ncol=2,nrow=0)
for(o in ls()[grepl('zzz',ls())]) config <- rbind(config, get(o))
write.csv(config,file=sprintf('%s/config.csv',rd_dir))


## Age bins, seperate config
ab_times <- data.frame(
                  tstart = c(0,29/30,6,12,24,36,48),
                  tstop  = c(29/30,6,12,24,36,48,60),
                  ab     = c('nn','pnn1','pnn2','1y','2y','3y','4y') )
ab_times$ab <- paste0(1:length(ab_times$ab),"_",ab_times$ab)
ab_times$tmid <- ab_times$tstart+(ab_times$tstop-ab_times$tstart)/2
write.csv(ab_times,file=sprintf('%s/ab_times.csv',rd_dir))



#######################################################
##### QSUB
rundates=run_date
rundates=c('2018_05_15_13_24_14',
           '2018_05_15_13_23_29',
           '2018_05_15_13_22_59',
           '2018_05_15_13_22_34')



# Loop through countries
countries <- read.csv(sprintf('%s/country_list.csv',root),stringsAsFactors=FALSE)$countries

for(run_date in rundates){
  for(country in countries){
    # make qsub
    qsub <- make_qsub(rtdir       = root,
                      rd          = run_date,
                      cores       = 20, #as.numeric(zzz_cores[1,2]),
                      memory      = 200,
                      proj        = 'proj_geospatial',
                      coderepo    = '/share/geospatial/royburst/code/sbh_cbh',
                      codescript  = 'CROSSVAL_holdout_run_parallel.R',
                      geo_nodes   = TRUE,
                      singularity = TRUE,
                      cntry       = country  )

    system(qsub)
  }
}
#
system('qstat')
#cd: substr(qsub,8,73)
