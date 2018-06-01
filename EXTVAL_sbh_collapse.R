## #######################################################
## AUTHOR: Roy Burstein
## DATE: MAY 2018, orig JANUARY 2018
## PURPOSE: This script combines adjusted SBH yearly data for each
#           SBH survey into one file and collapses onto the desired
#           years run_sbh_prep_models_predict.R.
#           Used to verify the method using external concurrent sbh data.
#           SBH - CBH paper version
#  OUPUT: A  prepped sbh combined and collapsed dataset
## #######################################################




## ##############################################################
## USER DEFINED GLOBALS FOR RUN DATE AND PERIOD TABULATION
user       <- Sys.info()['user']


# Define the sbh prediction run date to source individual prepped surveys from.
# This code looks for results in here: /share/geospatial/mbg/u5m/died/input/sbh/dated_sbh_model_prediction/
sbh_prep_rd <- "2018_05_15_14_40_13" 

# Start of each period.
period_start <- c(1990:2017)
period_label <- c(1990:2017)

## ##############################################################






## ##############################################################
## SETUP - LOOK FOR COMPLETED NID PREPS TO BRING IN, WARN ABOUT MISSING ONES

# libraries
library(data.table)
library(parallel)

# root dir for outputs
root <- '/share/geospatial/mbg/u5m/died/input/sbh/'

# load custom functions
source(sprintf('/share/code/geospatial/%s/u5m/bh_prep/sbh_utils.R', user))



# source a list of all the NIDs that should be in there:
nid_dirs <-  list.dirs(sprintf('%s/dated_sbh_model_prediction/%s',root,sbh_prep_rd))[-1]
all_nids <- gsub(sprintf('%s/dated_sbh_model_prediction/%s/',root,sbh_prep_rd),'',nid_dirs)


# look through each directory and search for a completed file, record where we have and dont have data
# TODO: Update this to look at if only some subsets finished
missing_prep_nids <- missing_prep_dirs <- as.character()
for(d in 1:length(nid_dirs)){
  if(!file.exists(sprintf('%s/fin.txt',nid_dirs[d]))) {
    missing_prep_nids <- c(missing_prep_nids,all_nids[d])
    missing_prep_dirs <- c(missing_prep_dirs,nid_dirs[d])
  }
}

nids_with_data <- all_nids[!all_nids %in% missing_prep_nids]
dirs_with_data <- nid_dirs[!nid_dirs %in% missing_prep_dirs]

if(length(missing_prep_nids) > 0) {
  message(sprintf('%i EXPECTED SBH PREP FILES ARE MISSING!/n',length(missing_prep_nids)))
  message(paste('THE FOLLOWING NIDS DO NOT HAVE PREP FILES:',paste(missing_prep_nids,collapse=', '),'/n'))
  message('WARNING: Please make sure all predict runs have finished,
          and if so and there are still missing files, check logs for errors.
          this code will continue to collapse and combine only those surveys
          which have prepped data for.')
}

## ##############################################################





## ##############################################################
## For those with data, collapse and summarize them

# write a function that loads data, collapses and summarizes
load_collapse_summarize <- function(dir){
  # store the nid
  nid <-  gsub(sprintf('%s/dated_sbh_model_prediction/%s/',root,sbh_prep_rd),'',dir)

  # load in the data, each subset file
  files <- list.files(dir, pattern = 'prepped_yearly_aggregation_draws_')
  d <- list()
  for(f in files)
    d[[f]] <- readRDS(sprintf('%s/%s',dir,f))
  d <- do.call('rbind', d)

  d[, period := findInterval(year_entering, period_start)]

  # keep only those with any estimates in a period greater than 0
  d <- subset(d, period > 0) # could also check based on svyyr

  # sometimes we only catch a few women at the last year because a couple ints were in jan
  # this leads to wonky estimates in the last year, so if the last year SS is tiny, drop it.
  # other than small ss, this is likely just some subset of the pop, so not representative
  tmp         <- d[,.(N=sum(EEB)),by=.(period)]
  outlier_val <- mean(tmp[period!=max(period)]$N) - 1.5 *sd(tmp[period!=max(period)]$N)
  if(tmp[period==max(period)]$N < outlier_val)
    d <- subset(d, period < max(d$period))

  if(nrow(d) == 0){ # if there are no data in the study period, end the function

    message(sprintf('FOR NID %s NO OBSERVATIONS WITHIN STUDY PERIOD. RETURNING BLANK DT',nid))
    return(data.table())

  } else { # if there are data within the correct period, then continue on with the function

    # remove year entering so can collapse on period instead
    d$year_entering <- NULL

    # collapse to the geo_unit-period-age bin: weighted mean each draw, sum EEB
    # NOTE, this is expecting draws to be of format VXXXX, V and up to 4 digits.
    collapse_on  <- colnames(d)[grep('V',colnames(d))] # columns for draw variables, (Must be V then some numbers)

    # removing geo_info so we simply collapse on country year
    d[ ,`:=`(geo_unit = NULL, point = NULL, longnum = NULL, latnum = NULL, location_code = NULL, shapefile = NULL)]

    # get columns to collapse on
    collapse_by  <- colnames(d)[grep('V|EEB|sum_wgt',colnames(d),invert=TRUE)]

    d <- d[, c(list(N = sum(EEB)), lapply(.SD, weighted.mean, sum_wgt)), .SDcols = collapse_on, by = collapse_by]


    # post-collapse 1qo ( 1 - prod(1-c(NN,PNN1,PNN2))) and 5q0 for each country
    collapse_ab_by <- collapse_by[collapse_by!='ab']

    # collapse the ages, use neonatal birth N as the sample size for the full period
    u5  <- d[, c(list(N = max(N)),
               lapply(.SD, function(x){1 - prod(1-x) } )),
               .SDcols = collapse_on, by = collapse_ab_by]
    inf <- subset(d, ab %in% c('1_NN','2_PNN1','3_PNN2'))[, c(list(N = max(N)),
                                                              lapply(.SD, function(x){1 - prod(1-x) } )),
                                                          .SDcols = collapse_on, by = collapse_ab_by]
    nn  <- subset(d, ab %in% c('1_NN'))
    nn[,  ab := 'NN']
    inf[, ab := '1q0']
    u5[,  ab := '5q0']
    d <- do.call('rbind', list(nn,inf,u5))

    # summarize draws, mean upper lower, 2.5% and 97.5%
    summ <- data.table(t(d[,apply(.SD,1,quantile,probs=c(0.5,0.025,0.975)),.SDcols=collapse_on]))
    colnames(summ) <- c('haz','haz_lower','haz_upper')

    # bind just the summaries back on
    d <- cbind(d[, c(collapse_by, 'N'), with = FALSE], summ)

    # return the dataset
    return(d)
  }
}

# apply over each survey and do the thing
message(sprintf('Loading, collapsing, and summarizing %i surveys.', length(dirs_with_data)))
sbh <- mclapply(dirs_with_data, load_collapse_summarize, mc.cores = 30)

# rbind over the list to get all the data into one frame
sbh <- do.call('rbind',sbh)
sbh$year <- period_label[sbh$period]


# median estimate of died for each row
sbh[, died := N*haz]

# clean up some variables to make sure the final prepped sbh dataset matches the final prepped cbh dataset
sbh[, data_type := 'sbh']
info <- fread('/share/geospatial/mbg/u5m/died/input/sbh/inputs/ihme_lc_info.csv')[,c('country','loc_id'),with=FALSE]
sbh <- merge(sbh, info, by='loc_id', all.x = TRUE)
sbh[,loc_id := NULL]


# print a quick summary of the data
print(sbh[, .(haz = weighted.mean(haz,N), N = sum(N)), by = .(ab,period)][order(ab,period)])
message(sprintf('There are %s approx. births in the SBH dataset', prettyNum(sum(sbh[ab=='NN']$N), big.mark = ',')))
message(sprintf('There are %s approx. deaths in the SBH dataset', prettyNum(sum(sbh$died), big.mark = ',')))
print(summary(sbh$haz_upper))

### ##############################################################





## ##############################################################
## SAVE

# Save a version
saveRDS(sbh, '/share/geospatial/royburst/sbh_cbh/sbh_verification/prepared_verification_sbh_tmp.RDS')
saveRDS(sbh, '/home/j/temp/royburst/sbh_cbh/verification_data/prepared_verification_sbh_tmp.RDS')

### ##############################################################
