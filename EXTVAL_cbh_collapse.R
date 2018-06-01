## ##############################################################
## Author: Roy Burstein
## Inputs: Child-Level CBH data prepared from Ubcov
## Outputs:Age/period-specific number entering bin and deaths as counts.
## Description: This script collapses CBH data User can input parameters for
##   period and age bin breakdowns. Survey weights will be used.
## Dependencies: data.table, Biograph, survival
## ##############################################################


## ##############################################################
## USER DEFINED GLOBALS FOR AGE BIN AND PERIOD TABULATION

## Age bins, define the start time (in months) of each age bin, and a label for each bin
ab_times <- data.frame( tstart = c(0, 1,  6, 12, 24, 36, 48),   # time start each age bin, in months
                        tend   = c(1, 6, 12, 24, 36, 48, 60),
                        ab     = c('NN','PNN1','PNN2','1yr','2yr','3yr','4yr') ) # age bin labels

## Time periods, define the starting year of each period
period_start <- c(1988:2017)
period_label <- c(1988:2017)

## set file sytem paths
# where to pull in the ubcov'd CBH data in from
dat_dir  <- '/home/j/WORK/11_geospatial/02_processed data/U5M/Global'
# where to ultimately save this
save_dir <- '/share/geospatial/mbg/u5m/died/input/data/'

## ##############




## ##############################################################
## SET UP
## Loading Packages
library(data.table)
library(Biograph)
library(survival)

## Loading in child-level CBH data
d <- setDT(readRDS(sprintf("%s/Global_CBH_1988_2017.RDS",dat_dir)))

# limiting to 1988 data
d <- subset(d, year >= 1988)
## ###############21331




## ##############################################################
## CLEAN DATA AND REPORT MISSINGNESS VALUES

## Standardize yearborn definition and check for missingness
d[, yrborn  := cmc_as_year(interview_date_cmc - birthtointerview_cmc)]
d[, yrborn2 := cmc_as_year(child_dob_cmc)]
d[is.na(yrborn), yrborn := yrborn2]
if(sum(is.na(d$yrborn)) > 0)
  message(sprintf('WARNING: There are %s rows (%f%%) with a missing yrborn variable.\n',
                  prettyNum(sum(is.na(d$yrborn)),big.mark=','),
                  sum(is.na(d$yrborn))/nrow(d)*100))


## Standardize year of interview definition and check for missingness
d[, yrint := cmc_as_year(interview_date_cmc)]
if(sum(is.na(d$yrint)) > 0)
  message(sprintf('WARNING: There are %s rows (%f%%) with a missing yrint variable.\n',
                  prettyNum(sum(is.na(d$yrint)),big.mark=','),
                  sum(is.na(d$yrint))/nrow(d)*100))


## Standardize aod definition and check for missingness and oddities
## Also, define variable 'alive', will get used as a censoring indicator
setnames(d,old='child_age_at_death_months',new='aod')
# NOTE FOR LATER: if want ENN will need to use raw AOD to get daily and weekly precision, for now just count NN as 1 month
if(as.numeric(names(sort(-table(d$aod))))[1] != 6000)
  stop('The mode of AOD was not 6000, check that "alive" was coded as 6000 as expected')
d[,alive := aod == 6000]  # if alive at time of survey aod is coded as 6000, so consider them alive
d[,alive := aod >= 60]   # Consider all those who died after 60 months as 'alive' for our purposes of studying under-5 mortality. Equivalent to being right-censored from this perspective
d$alive[is.na(d$child_alive)] <- NA # if child_alive is unknown in the ubcov data, keep as NA
d[,died  := alive == FALSE] # if not alive, then they have died
if(sum(is.na(d$aod))>0) message(sprintf('WARNING: A total of %s rows (%f%%) have missing age of death (aod).\n',
                                        prettyNum(sum(is.na(d$aod))),
                                        sum(is.na(d$aod))/nrow(d)*100))
if(sum(is.na(d$alive))>0) message(sprintf('WARNING: A total of %s rows (%f%%) have missing alive/died indicators.\n',
                                          prettyNum(sum(is.na(d$alive))),
                                          sum(is.na(d$alive))/nrow(d)*100))
if(sum(d$died,na.rm=TRUE)==0)
  stop('There are no dead children in this dataset. Very unlikely. Please confirm.')


## Standardize child age definition and check for missingness (record age in months)
#  These time-specific indications of alive and dead are used later on the survival:: reshape.
d[,childage := round((yrint - yrborn)*12)]         # if alive, record their age
d$childage[d$died==TRUE & !is.na(d$died)] <- d$aod[d$died==TRUE & !is.na(d$died)]    # if died, record their age at death
d[,childage := childage + 0.0001]      # add a tiny bit to help with discrete age category definitions
if(any(d$childage<0)) message('WARNING: There are negative aged children in the data! Something could be wrong with dates.')
if(min(d$childage,na.rm=TRUE)>0.0001) message(sprintf('WARNING: The youngest child in this data is %f, this number should be 0.',min(d$childage)))
if(max(d$childage,na.rm=TRUE)>42*12) message('WARNING: The oldest child is >42yo, this is highly unlikely as max mother age responding is 49. Confirm.')


## if an entire survey is missing weights, then we give them a 1 weight.
#   This follows on the explicit choice to include surveys without weights as unweighted
#   The only non-weighted rows remaining after this will be polygon rows in surveys with
#   weights, theyll be dropped later on
d[nid==26661, weight := 1] # for this eth GF survey for some reason the weights are all 0..
nowgt <- d[,.(pct_unweighted = sum(is.na(weight))/.N), by = nid]
nowgt <- subset(nowgt, pct_unweighted == 1)$nid
if(length(nowgt) > 0) {
  message(sprintf('WARNING: %i surveys have no weight info and will be used as unweighted.',length(nowgt)))
  d$weight[d$nid %in% nowgt] <- 1
}

## Report how many children are being dropped and then drop them before continuing
d[, must_drop := is.na(weight) | is.na(childage) | childage > 41*12 |
    childage < 0 | is.na(aod) | is.na(alive) | is.na(yrborn)]
if(sum(d$must_drop) > 0)
  message(sprintf('WARNING: Dropping %s rows (%f%%) from the final dataset due to missingess in these variables.',
                  prettyNum(sum(d$must_drop)),sum(d$must_drop)/nrow(d)*100))
d <- subset(d, must_drop == FALSE)


## Make the DT a bit smaller to speed up the reshape coming up next
# keep only the variables we need
d <- d[, c('nid','source','country','year','strata','weight',
           'yrborn','aod','childage','died','yrint'),
       with = FALSE]

# keep only children born within the range of dates we want
d <- subset(d, yrborn > (min(period_start) - 6))

message(sprintf('INFORMATION ON %s CHILDREN FROM %i SURVEYS TO BE USED.',
                prettyNum(nrow(d),big.mark=','),length(unique(d$nid))))

## ################





## ##############################################################
## RESHAPE AND IDENTIFY AGE BINS AND PERIODS

## make the data long, such that each row is a child-agebin (Age bins are defined at the top of this script)
cbh <- setDT(survival::survSplit(Surv(childage, died) ~., cut = ab_times$tstart, data=d))
if(sum(d$died) != sum(cbh$died))
  message('The number of children died pre and post data reshpae (age-wise lengthing) are not the same. Investigate.')


# identify the year each discrete age bin was entered by each child
cbh <- merge(cbh, ab_times, by = 'tstart') # get age bin names
cbh[, year_entering := yrborn + tstart/12 ]

# check how many are observed alive.. (censored) and are just being included in the probability demoninator
nrow(cbh[(childage < tend & died == 0)])/nrow(cbh) # 4% of our obs.. but also would need to delete deaths of children who would be censored if they had lived to the censoring point:
nrow(cbh[(round((yrint-year_entering)*12+tstart) < tend & died == 1)])/nrow(cbh)
# all possibly censored on slate to be dropped under this framework:
cbh[,censdrop:=(round((yrint-year_entering)*12+tstart) < tend )]

# identify the period bin we wish to use (Period bins are defined at the top of this script)
cbh[, period := findInterval(year_entering, period_start)]
# drop any info on age bins for time period 0, which is before the first period
cbh <- subset(cbh, period > 0)

# associate a year with a period
setnames(cbh, 'year', 'svyyr')
cbh$year <- period_label[cbh$period]


# print out the mean hazards by age bin and period for this survey, add this to report later
print(cbh[, .(hazard_of_death = sum(died)/.N, SS = .N), by = .(ab,year)])
print(cbh[censdrop==FALSE, .(hazard_of_death = sum(died)/.N, SS = .N), by = .(ab,year)])

## ################





## ##############################################################
## COLLAPSE, this part replaces seegMBG::periodTabulate()

##  collapse to urvey year, use survey weighted mean (make sure weights are normalized to 1000 per survey)
cbh_final <- cbh[censdrop==FALSE, .(haz = weighted.mean(died, w = weight), N = .N),
                 by = .(country,nid,source,svyyr,period,year,ab)]
cbh_final[, died := haz*N]
cbh_final[, data_type := 'cbh']

# subset these to the 15 years prior to the survey
cbh_final <- subset(cbh_final, svyyr - year < 16)

# sometimes we only catch a few women at the last year because a couple ints were in jan
# this leads to wonky estimates in the last year (ie year of survey. well drop these, so the first est comes from one year before survey)
cbh_final[, mxy := year == max(year), by = nid]
cbh_final <- subset(cbh_final, mxy == FALSE)
cbh_final[, mxy := NULL]

## ########

### ADD YEAR VARIABLE


## ##############################################################
## SAVE THE PREPARED CBH DATASET

saveRDS(cbh_final, '/home/j/temp/royburst/sbh_cbh/verification_data/prepared_verification_cbh_1988_2017_censfix.RDS')

## ########
