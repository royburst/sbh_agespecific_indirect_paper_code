message(' ... ... ... ... ... ... ... ... ... ... ')
message(' ... ... BEGIN VALIDATION RUN ... ... ...')
message('Run notes:')
print(R.Version())
print(Sys.info())
message(' ... ... ... ... ... ... ... ... ... ... ')
message(' ... ... ... ... ... ... ... ... ... ... ')

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~ SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Load in options from qsub call
run_date               <- as.character(commandArgs()[4])
cntry                  <- as.character(commandArgs()[5])

message(sprintf('\nRUN_DATE: %s',run_date))
message(sprintf('\nCOUNTRY: %s',cntry))




######################################################################
message('\nSETUP ... ')

## Get output directory
root   <- '/share/geospatial/royburst/sbh_cbh/output'

## Load libraries (should be running from singularity image!)
packlist <- c('data.table','ggplot2','parallel','tictoc','raster','rgeos',
              'seegMBG','survival','rms','flexsurv','magrittr','scales','sp','maptools',
              'gridExtra','lme4','grid','INLA','Biograph')
for(pack in packlist) {
  message(paste0('Loading from user list: ', pack))
  library(pack,character.only=TRUE)
}

# ill need some MBG things
repo <- sprintf('/share/code/geospatial/%s/mbg/','royburst')
setwd(repo)
for(funk in list.files(recursive=TRUE,pattern='functions')){
  if(length(grep('central',funk))!=0){
    message(funk)
    source(funk)
  }
}

# load custom functions
setwd('/share/geospatial/royburst/code/sbh_cbh')
source('utils.R')
message(' ... Custom functions sucessfully loaded')

#
make_output_dirs(root = root, run_date = run_date, country = cntry)
rd_dir  <- sprintf('%s/%s', root, run_date)
rdc_dir <- sprintf('%s/%s', rd_dir, cntry)
dat_dir  <- '/home/j/WORK/11_geospatial/02_processed data/U5M/Global'

## load configs, which were saved in launch script
# ab times
ab_times <- fread(sprintf('%s/ab_times.csv',rd_dir))

# config values into environment
load_config(rd_dir)




######################################################################
## LOAD IN THE DATA
message('\nLOAD IN CBH-SBH COMBINED DATASET ... ')
tic.clearlog()
tic('data_loading')

data <- readRDS("/home/j/WORK/11_geospatial/02_processed data/U5M/Global/CBH_with_SBH_1988_2017.RDS")
data <- data.table(data)
data <- merge_ihme_loc(data)
data <- subset(data, survey %in% c('MACRO_DHS','MACRO_MIS','NIC/DHS_ENDESA'))
data <- subset(data,! nid %in% 20786) # sen out
data$nid <- as.factor(data$nid)
message(sprintf(' ... Full database with %i rows sucessfully loaded ',nrow(data)))
toc(log = TRUE) # end data_loading tic




######################################################################
# LOAD IN SDI INFO
message('\nLOAD SDI INFORMATION ... ')

sdi <- fread('/share/geospatial/royburst/sbh_cbh/data/sdi.csv')[,c('location_id','year_id','sdi'),with=FALSE]
setnames(sdi,c('location_id','year_id'),c('loc_id','year_entering'))
sdi[, year_entering := year_entering-base_year]
message(' ... SDI information sucessfully loaded')




######################################################################
message('\nSORT OUT NIDS FOR TEST AND TRAIN')
nids <- unique(data[,c('nid','year','country','survey','loc_id','MAP_gbdregion'),with=FALSE])
nids$nid <- as.factor(nids$nid)

# get training NIDS based on config extent
if(training_extent == 'global'){
  usenids <- nids$nid
} else if(training_extent == 'region'){
  region <- nids$MAP_gbdregion[nids$country == cntry][1]
  usenids <- nids[MAP_gbdregion == region]$nid
} else if(training_extent == 'country'){
  usenids <- nids$nid[nids$country == cntry]
  # if country is te, do not allow for country or country_ab random effects in the formula
  formula <- gsub('\\+ s\\(country_ab, bs = \'re\'\\)','',formula)
  formula <- gsub('\\+ s\\(country, bs = \'re\'\\)','',formula)
} else if(training_extent == 'survey'){
  stop('training_extent == survey is not yet supported')
}

# get test NID based on strategy, only last_year in country currently supported ho strategy
if(ho_strat == 'last_year'){
  nid_test <- nids[nids$country == cntry][order(-year)]$nid[1]
  nid_train <- usenids[usenids != nid_test] # remove test from train list
} else {
  stop('last_year is the only currently supported ho_strat')
}
message(sprintf(' ... Test nid: %s', nid_test))

# write out the test train NID identifiers, in case we need them later
nids$train <- NA
nids$train[nids$nid %in% nid_train] <- TRUE
nids$train[nids$nid %in% nid_test ] <- FALSE
nids <- subset(nids, !is.na(train))
write.csv(nids, file = (sprintf('%s/nids.csv',rdc_dir)))




######################################################################
# Subset out the data
message(' \nDATA SUBSET TO USENIDS ONLY ....')
tic('data_prep')
d <- subset(data, nid %in% usenids)

# clean up data
message('Cleaning up data .... \n')
d <- clean_dhs(d)
d <- subset(d,childage>=0)
d$yrborn_orig <- d$yrborn




######################################################################
message('\nCOVARIATES TO UNIT SCALE ....')
# parse out all the fixed effects variables from the gam formula
vars <- get_fe_from_gamformula(formula)

# centrescale all covariates and save for later
cs   <- get_cs_table(d,vars, skipvars = c('sdi','birthorder'))
for(v in cs$name) {
  d[[paste0(v,'_orig')]] <- d[[v]]
  if(v %in% colnames(d)) # skips SDI here since its not loaded in yet
    d[[v]] <- centscale(d[[v]],v,reverse=FALSE)
}




######################################################################
message('\nRESHAPE LONG AND FORMAT FOR DTSA ALONG SPECIFIED AGE BINS .... ')
# use survival package handy tools for this
dl <- data.table(survSplit(Surv(childage, died) ~.,cut=ab_times$tstart,data=d))
dl <- merge(dl,ab_times,by='tstart') # add on age bin information
dl[, year_entering := round(yrborn_orig + tstart/12,0)] # identify the year entering each bin

# drop observations that are censored within the most recent bin
# time entering + bin length must be less than survey year
dl[, cens := (yrborn_orig + tstart/12) + (tstop - tstart)/12 > (cmc_as_year(interview_date_cmc) - base_year) ]
message(sprintf('%i are censored (%fpct)',sum(dl$cens,na.rm=TRUE),round(sum(dl$cens,na.rm=TRUE)/nrow(dl)*100,2)))
dl <- subset(dl, cens == FALSE)

dl <- merge(dl,sdi,by=c('loc_id','year_entering'),all.x=TRUE) # add SDI
dl[,dum := 1] # make a dummy which we use later in period aggregation steps
dl[,fyrint := floor(yrint),by=nid] # floor the year of interview


# anything used as a random effect should be labeled as a factor
dl$mother_id <- as.factor(dl$mother_id)
dl[, country_ab := as.factor(paste0(country,'_',ab))]
dl$country_orig <- dl$country
dl$country   <- as.factor(dl$country)
dl$nid_yrborn   <- as.factor(paste0(dl$nid,'_',round(dl$yrborn,2)))


# backfill pre 1970, not used
sdib <- dl[,.(sdib = min(sdi,na.rm=TRUE)),by=loc_id]
dl <- merge(dl,sdib,by='loc_id')
dl$sdi[is.na(dl$sdi)] <- dl$sdib[is.na(dl$sdi)]




######################################################################
message('\nSPLIT OUT DATASET AS TRAINING AND TESTING COMPONENTS .... ')
# split training and testing datasets
if('oos' %in% eval(parse(text=predict_for))){
  message(' ... Splitting dl into test and train')
  dl_test   <- dl[nid %in% nid_test]
  dl_train  <- dl[nid %in% nid_train]
} else {
  message(' ... Out of sample toggle was off so not splitting data')
  dl_test  <- dl_train <- dl
}
toc(log = TRUE) # end data_prep tic

# scale weights so they sum to the same sample size, because of the way mgcv uses weights
# basically the total N should still be the same so we dont change the total contribution to the likelihood
dl_train[, weight := weight/sum(weight,na.rm=TRUE)*(.N), by = nid]




######################################################################
message(sprintf('\nRUNNING GAM (mgcv::bam) MODEL ON %i CORES AND SAVING OBJECT... ',cores))
# build cluster with X cores and run a bam model
tic('model_fit')
cl  <- makeCluster(cores)
mod <- mgcv::bam(as.formula(formula), data = dl_train, family = 'binomial', cluster = cl , weight = weight)
stopCluster(cl)
toc(log = TRUE) # end model_fit tic
saveRDS(mod, file = sprintf('%s/model_object.RDS',rdc_dir))

print(summary(mod))





######################################################################
message('\nPROCEEDING TO PREDICTION ... ')
message(sprintf(' ... predict_for: %s',predict_for))
tic('prediction')

# predict now, either in or out of sample
#for(samp in eval(parse(text=predict_for))) {
samp <- 'oos'
message(paste0(' ... Predicting for, ',samp,' ... \n\n'))
oos <- ifelse(samp == 'oos', TRUE, FALSE)

# keep either in or out of sample data
if(samp=='oos') nid_keep <- nid_test   else  nid_keep <- nid_train

# data to predict  subsetted to kept train or test nids (depending on samp)
data_to_pred <- d[nid %in% nid_keep]

# make new data
message(' ... ... Making new (prediction frame) data, at mother level ')
nd  <- get_sbh_mother_info(data_to_pred, variables = vars, form = NULL)

# need mothers age at int to expand the frame porperly
nd$mothage_atint_orig <- centscale(var=nd$mothage_atint,varname='mothage_atint',reverse=TRUE)

message(' ... ... Expanding the frame to be by child-age')
nd  <- expand_prediction_frame(nd, mothageatintvar = 'mothage_atint_orig', abt = ab_times) #_orig')
message(sprintf(' ... ... ... Prediction frame with %i rows created',nrow(nd)))

# add SDI
message(' ... ... Adding year-varying covariates to child-age prediction frame')
nd <- merge(nd,nids[,c('nid','loc_id')],by='nid',all.x=TRUE)
nd[, year_entering := floor(yrborn + tstart/12)]

message(' ... ... ... SDI')
nd <- merge(nd,sdi,by=c('year_entering','loc_id'),all.x=TRUE)
nd <- subset(nd,year_entering < (2017-base_year)) # limit of SDI data
sdib <- nd[,.(sdib = min(sdi,na.rm=TRUE)),by=loc_id]   # backfill pre base_year
nd <- merge(nd,sdib,by='loc_id')
nd$sdi[is.na(nd$sdi)] <- nd$sdib[is.na(nd$sdi)]


# Re-scale all variables
message(' ... ... Center scaling imputed variables in new prediction frame')
for(v in c('yrborn','mothage_atbirth')) {
  nd[[paste0(v,'_orig')]] <- nd[[v]]
  nd[[v]] <- centscale(nd[[v]],v)
}

# make an orig ceb variable for naming convention, used later in weighting
if(!'ceb_orig' %in% colnames(nd)) nd$ceb_orig <- nd$ceb

# add on pct_ceb, based on MAP distributions
message(' ... ... ... Also imputing pct_ceb and birthorder variables')
map    <- load_MAP_distributions(gbdregions=data_to_pred$MAP_gbdregion[1])
nd     <- get_pct_ceb()

# centrescale pct_ceb if needed
nd$pct_ceb <- centscale(var=nd$pct_ceb_orig,varname='pct_ceb')
nd$birthorder_orig <- nd$pct_ceb_orig*nd$ceb_orig
nd$birthorder <- centscale(var=nd$birthorder_orig,varname='birthorder')

# DO A CHECK THAT ALL THE CENTERSCALING WORKED.
message(' ... ... ... Just a check to see that model data and prediction frame vars in the same range, what with all the center scaling')
for(v in vars){
  message(paste0('\n',v))
  message('Model Frame:'); print(summary(mod$model[[v]]))
  message('Prediction Frame:'); print(summary(nd[[v]]))
}

# make sure we also have country ab
nd[, country_ab := as.factor(paste0(country,"_",ab))]
nd$nid_yrborn   <- as.factor(paste0(nd$nid,'_',round(nd$yrborn,2)))

# rename weight survey weight to avoid confusion with EEB later
setnames(nd, 'weight', 'svywgt')

# !!!! TEST ZONE
if(TRUE == FALSE){
  for(i in 1:10) message(' !!!!!!!! WARNING. YOU ARE IN TEST MODE. REMOVING ALL BUT 100 MOTHERS !!!!!!!! ')
  ndo <- copy(nd)
  nd <- nd[mother_id %in% sample(unique(nd$mother_id),100)]
}
message(' ... ... Completed setting up prediction frame')




##################################################################
message(' ... ... Making individual level predictions on newdata')

res <- predict_Dsurv(model     = mod,
                     new_data  = nd,
                     yrbornvar = 'yrborn_orig',
                     ndraws    = ndraws,
                     summarize = FALSE,
                     ncores    = cores)

message(' ... ... Prediction is complete, aggregating period estimates')
# agregate over period using MAP to get exp number entering each bin
aggresl <- aggregate_period_estimates(est           = res,
                                     map             = map,
                                     yrbornvar       = 'yrborn_orig',
                                     cebvar          = 'ceb_orig',
                                     mothageatintvar = 'mothage_atint_orig',
                                     TESTNOMAP       = FALSE,
                                     use_svy_weight  = TRUE,
                                     yrs_to_agg      = yrs_to_agg)

# check to see if we have draw level aggregates
if(class(aggresl) == 'list'){
  draws_agg <- aggresl$draws
  draws_agg[, year := year_entering+base_year]; draws_agg[,year_entering := NULL]
  aggres    <- aggresl$res
  drawlevel <- TRUE
} else {
  aggres <- aggresl
  drawlevel <- FALSE
  draws_agg <- 0
}

# set up in normal year spaces
aggres[, year := year_entering+base_year]; aggres[,year_entering := NULL]

message(' ... ... Also, tabulate raw hazards from training and testing data for comparison:')
# collapse to get empirical hazard estimates from actual data for comparison
# collapse on age bin and year entering the bin to get an empirical period estimate for that bin
aggdat <- tabulate_raw(earliest_year  = 1990,
                       dl_train       = dl_train,
                       dl_test        = dl_test,
                       oos            = oos,
                       yrs_to_agg     = yrs_to_agg,
                       use_svy_weight = TRUE)
aggdat[,year_entering := NULL]


message(' ... ... Bind together raw tabulates and estimates to compare')
# make a comparison of period agg datasets
comp <- na.omit(merge(aggres,aggdat,by=c('ab','year')))
setnames(comp,c('haz.x','haz.y'),c('haz_est','haz_dat'))


message(' ... ... Getting 5q0 estimates out as well for comparison')
## get q5 using draws, if they are available
q5 <- get_q5(comp, c('haz_est','haz_lower','haz_upper','haz_dat'))
if(drawlevel) {
  draw_cols <- colnames(draws_agg)[!colnames(draws_agg) %in% c('ab','year','train')]
  draws_agg$train <- FALSE
  q5draw <- get_q5(draws_agg, draw_cols)
  q5draw <- cbind(q5draw[,c('ab','year','train'),with=FALSE],
                  data.table(t(q5draw[,apply(.SD,1,quantile,probs=c(0.5,0.025,0.975)),.SDcols=draw_cols])))
  colnames(q5draw) <- c('ab','year','train','haz_estd','haz_lowerd','haz_upperd')
  q5 <- merge(q5, q5draw, by=c('ab','year','train'),all.x=TRUE)
  for(v in c('haz_est','haz_lower','haz_upper' )){
    q5[[v]][!is.na(q5[[paste0(v,'d')]])] <- q5[[paste0(v,'d')]][!is.na(q5[[paste0(v,'d')]])]
    q5[[paste0(v,'d')]] <- NULL
  }
}


if(oos){
  message(' ... ... Calculating some predictive validity metrics')
  # get some metrics of predictive validity and save
  pv <- get_oos_pv_summary_metrics(rbind(comp,q5))
  write.csv(pv, file = sprintf('%s/oos_pv.csv',rdc_dir))
}


message(' ... ... Saving output datasets')
res <- res[['new_data']]
for(dat in c('q5','comp','aggres','aggdat','draws_agg','res','cs'))
  write.csv(get(dat), file = sprintf('%s/%s.csv',rdc_dir,dat))


message(' ... ... Making plots')
dashplots(filename=sprintf('%s/dashplot2.pdf',rdc_dir),oos=oos,form = as.formula(formula))


## Q5 PLOTS
# load in brass estimates
brass <- fread('/share/geospatial/royburst/sbh_cbh/data/brass_q5_estimates.csv')
brass <- brass[nid %in% as.integer(as.character(nid_keep)),][,c('q5','time'),with=FALSE]


# get IHME q5
library(readstata13)
ihm   <- setDT(read.dta13('/home/j/WORK/02_mortality/02_inputs/03_5q0/sbh/FINAL_ESTIMATES/EST_DHS_v5Q0_IHME.dta'))
nidyr <- nids$year[nids$train==FALSE]
ihm <- subset(ihm, iso3 == cntry & source_date >= (nidyr-1) & source_date <= (nidyr+1)) # THIS IS JENKY! but they dont use nids!
setnames(ihm,'t','time')
ihm$nid <- nids$nid[nids$train==FALSE]
ihm[,q5:=q5/1000]
write.csv(ihm, file = sprintf('%s/ihme_estimates.csv',rdc_dir))


# plot q5 info
if(oos) q5sub <- q5[train==FALSE,]  else q5sub <- q5[train==TRUE,]
q5plot(q5dat = q5sub,filename = sprintf('%s/q5plot.pdf',rdc_dir), brassdat = brass, ihmedat = ihm)

#}






########################################################################
# saving individual, cluster, ad1, ad2 level comparisons of IS and OOS
message('... ... Saving cluster level comparisons of IS and OOS data') # from ./save_individual_level_comparison.R
keep_nid <- nids$nid[nids$country == cntry]    #==FALSE] # identify the oos dataset
raw <- data[nid %in% keep_nid]
latlon <- setDT(subset(readRDS(sprintf("%s/Global_CBH_1988_2017.RDS",dat_dir)), nid %in% keep_nid))
latlon <- unique(latlon[,c('cluster_number','latnum','longnum','location_code','shapefile','nid'), with=FALSE])
latlon[,cluster_number := as.numeric(as.character(cluster_number))]
latlon[,nid := as.character(nid)]
raw[,nid := as.character(nid)]

#### prep the CBH data
## collapse raw cbh to the cluster level
rawi <- clean_dhs(raw)
rawi <- merge_ihme_loc(rawi)
rawi[,yrborn_orig := yrborn]
rawi <- subset(rawi,childage>=0)

# sort out point polygon ids
rawi <- merge(rawi,latlon,by=c('cluster_number','nid'),all.x=TRUE)
rawi[,point:=!is.na(latnum)]
rawi[point==FALSE, cluster_no:=paste0(location_code,shapefile)]
rawi[point==TRUE , cluster_no:=paste0(cluster_number)]

# make data long by agebin
rawi <- data.table(survSplit(Surv(childage, died) ~., cut=ab_times$tstart, data=rawi))

# bin by 5-year entering the age bin
rawi <- merge(rawi,ab_times,by='tstart')
rawi[, year_entering := round(yrborn_orig + tstart/12,0) + base_year]
#rawi[,period:=findInterval(year_entering, yrbins)] # cmt out to keep at year resolution for now
rawi <- subset(rawi, yrint - yrborn <= 24) # make sure no births prior to 24 years are in here from CBH to match SBH processing
#rawi <- subset(rawi, period > 0)           # remove anything pre earliest year

# drop observations that are censored within the most recent bin
# time entering + bin length must be less than survey year
rawi[, cens := (yrborn_orig + tstart/12) + (tstop - tstart)/12 > (cmc_as_year(interview_date_cmc) - base_year) ]
message(sprintf('%i are censored (%fpct)',sum(rawi$cens,na.rm=TRUE),round(sum(rawi$cens,na.rm=TRUE)/nrow(dl)*100,2)))
rawi <- subset(rawi, cens == FALSE)


# collapse: tabulate the cbh (unweighted for now)
cbh <- rawi[, .(entering = .N, died = sum(died)), by = .(cluster_no,year_entering,ab,nid,longnum,latnum,shapefile,location_code)]
#cbh <- subset(cbh, !is.na(cluster_no))


#### prep the SBH data
# do some archaeology based on mother id to get cluster_number
res[,V1:=NULL]
resi <- merge(res, unique(rawi[,c('mother_id','cluster_number'),with=FALSE]),
              by='mother_id',all.x=TRUE)
#resi <- subset(resi, !is.na(cluster_number)) # WORRIED SOME OF THESE ARE GETTING DROPPED UNNCESSARILY!

resi <- merge(resi,latlon,by=c('cluster_number','nid'),all.x=TRUE)
resi[,point:=!is.na(latnum)]
resi[point==FALSE, cluster_no:=paste0(location_code,shapefile)]
resi[point==TRUE , cluster_no:=paste0(cluster_number)]



# get number entering and died at each row
resi    <- aggregate_period_estimates(est               = resi,
                                      map               = map,
                                      yrbornvar         = 'yrborn_orig',
                                      cebvar            = 'ceb_orig',
                                      mothageatintvar   = 'mothage_atint_orig',
                                      TESTNOMAP         = FALSE,
                                      yrs_to_agg        = 1,
                                      return_individual = TRUE)
#resi <- resi[,period:=findInterval(year_entering+1970, yrbins)]
#resi <- subset(resi, period > 0) # remove anything pre earliest year
resi[, year_entering := year_entering + base_year]
sbh <- resi[, .(entering = sum(EEB), died = sum(haz*EEB)), by = .(cluster_no,year_entering,ab,nid,longnum,latnum,shapefile,location_code)]

### save prepped data
cbh$source <- 'raw cbh'
sbh$source <- 'indirect sbh'
resfinal <- rbind(cbh,sbh)
latlon[,nid := as.character(nid)]
resfinal[,N:=entering]
setnames(resfinal,c('longnum','latnum'),c('longitude','latitude'))
resfinal[, age_bin := as.numeric(substr(ab,1,1))]
resfinal$yrbinstart <- resfinal$year_entering #yrbins[resfinal$period]
resfinal <- merge(resfinal,nids[nid%in%keep_nid][,c('nid','train'),with=FALSE],by='nid')
resfinal[,latitude := as.numeric(latitude)]
resfinal[,longitude := as.numeric(longitude)]
# poly resample - Keeps breaking, ignore it for now

write.csv(resfinal,sprintf('%s/cbh_sbh_compare_cluster_level.csv',rdc_dir))



###################################
### Save admin 1 level as well
message('... ... Saving Ad1 comparisons of IS and OOS data') # from ./save_individual_level_comparison.R
a1 <- readRDS('/share/geospatial/rds_shapefiles/admin2013_1.rds'); a1 <- a1[a1@data$COUNTRY_ID==cntry,]

# since poly resample didnt work, just get centroid of it and use that to id adm1
for(shp in na.omit(unique(resfinal$shapefile))){
  if(shp!=''){
    s <- readRDS(sprintf('/share/geospatial/rds_shapefiles/%s.rds',shp))
    sc <- setDT(data.frame(cbind(coordinates(gCentroid(s,byid=TRUE)),location_code=s@data$GAUL_CODE,shapefile=shp)))
    resfinal <- merge(resfinal,sc,by=c('shapefile','location_code'),all.x=TRUE)
    resfinal$longitude[!is.na(resfinal$x)] <- resfinal$x[!is.na(resfinal$x)]
    resfinal$latitude[!is.na(resfinal$y)]  <- resfinal$y[!is.na(resfinal$y)]
    resfinal[,x:=NULL]
    resfinal[,y:=NULL]
  }
}
resfinal <- subset(resfinal, !is.na(latitude) | !is.na(longitude))
pts <- SpatialPoints(cbind(resfinal$longitude, resfinal$latitude),
                     proj4string = CRS(proj4string(a1)))
adagg <- cbind(resfinal,over(pts,a1))
adagg <- adagg[,.(died=sum(died),N=sum(N)),
               by=.(nid,train,source,age_bin,ab,year_entering,GAUL_CODE,NAME)]
adagg[, mr := died/N]

write.csv(adagg,sprintf('%s/cbh_sbh_compare_ad1_level.csv',rdc_dir))





###################################
### collapse ad1 to 5 year bins and plot IS vs OS survey only
message('... ... Plotting AD1 comparisons') # from ./save_individual_level_comparison.R
yrbins <- seq(1993,2013,by=5)
ap <- adagg[nid == nids$nid[nids$train==FALSE]]
ap <- ap[,period:=findInterval(year_entering, yrbins)]
ap <- subset(ap, period > 0)
ap$yrbinstart <- yrbins[ap$period]
ap <- ap[,.(died=sum(died),N=sum(N)),
            by=.(nid,train,source,age_bin,ab,yrbinstart,GAUL_CODE,NAME)]
ap[, mr := died/N]

# scatter
apt <- merge(ap[source=='indirect sbh'],ap[source=='raw cbh'],by=c('ab','yrbinstart','GAUL_CODE'))
pdf(sprintf('/%s/test.pdf',rdc_dir),height=12,width=8)
    ggplot(apt,aes(x=mr.x,y=mr.y,col=factor(yrbinstart),size=N.y)) + theme_bw() + ylab('CBH') + xlab('OOS SBH')+
             geom_point()+geom_abline(intercept=0,slope=1,col='red') + facet_wrap(~ab,ncol=1) +
              ggtitle(paste('Admin 1 level estimates comparison,',cntry,', NID:',nids$nid[nids$train==FALSE]))
#dev.off()

# place these back on maps to plot
colz8 <- c('#2A161E','#312B3B','#264451','#1F5E56','#41744B','#79843C','#BB8C43','#F88F6E')
colz <- scales::seq_gradient_pal(colz8)(seq(0,1,length.out=100))
#pdf(sprintf('/%s/cbh_sbh_ad1_compare.pdf',rdc_dir),
#    height=length(yrbins)*4,width=8)
par(mfrow=c(length(yrbins),2))
for(a in sort(unique(ap$ab))){
  tmp1 <- subset(ap, ab == a)
  brks <- cut(tmp1$mr,breaks = 8)
  leglab <- c('0',rep('',6),as.character(round(max(tmp1$mr),3)))
  tmp1$col <- colz8[as.numeric(cut(tmp1$mr,breaks = 8))]
  for(y in 1:length(yrbins)){
    for(so in unique(ap$source)){
      tmp <- subset(tmp1, yrbinstart == yrbins[y] & ab == a & source == so)
      p <- merge(a1,tmp,by='GAUL_CODE')
      plot(p, col = p@data$col, main = paste(a,yrbins[y],so))
      if(y==1) legend("bottom", fill = colz8, legend = leglab, col = colz8, box.lwd = 0)
    }
  }
}
dev.off()




################################################################################################
toc(log = TRUE) # end prediction tic






######################################
message('SAVING TIC TOC LOG .... \n')
# grab toc log and format it nicely and save
tl <- t(matrix(unlist(strsplit(tl,': ')),nrow=2))
write.csv(tl, file = sprintf('%s/timing_log.csv',rdc_dir))



######################################
message('SAYING GOODBYE .... \n')
#
