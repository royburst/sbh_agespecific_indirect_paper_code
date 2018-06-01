## Roy Burstein
## Functions used in the analysis






# Work with time as decimal years (with 0 = Jan 1970)
# Helper functions to convert from cmc (century month code)
cmc_yr_convert <- function(x, zero_year = 1970) {
  require(Biograph)
  cmc_as_year(x) - zero_year
}



# helper find mode
get_mode <- function(x,top=1) as.numeric(names(sort(-table(x))))[1:top]



# helper trim whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)



# helper grab first observation
first <- function(x) x[1]



# helper get J root depending on which OS you are using
j <- function() ifelse(Sys.info()[['sysname']]=='Windows','J:/','/home/j/')



# Write notes to myself for each model run date
write_rd_notes <- function(note, rd = run_date, rt = root){
  run_notes <- read.csv(sprintf('%s/run_date_notes.csv',rt), stringsAsFactors = FALSE)[,2:3]
  dup_rd <- ''
  if(rd %in% run_notes$run_date){
    while(! dup_rd %in% c('Y','N'))
    dup_rd <- readline(sprintf('run_date %s is already in the run_notes file. Wish to overwrite? (Y/N)', rd))
  }
  if(dup_rd == 'N') {
    stop('Will not overwrite dup, maybe grab a new run_date?')
  } else if (dup_rd == 'Y') {
    run_notes <- run_notes[run_notes$run_date != rd, ]
  }
  run_notes <- rbind(cbind('run_date'=rd,'note'=note),run_notes)
  message('Writing new note.')
  write.csv(run_notes,sprintf('%s/archive_run_date_notes/run_date_notes_%s.csv',rt,rd))
  write.csv(run_notes,sprintf('%s/run_date_notes.csv',rt))
}






# format some things in DHS data for a smooth analysis later
clean_dhs <- function(data,
                      zero_year  = 1970, # year we count as zero
                      censor_age = 60,   # age at which we censor observations (considered outside the study period)
                      alive_code = 6000  # numeric age code for alive at time of survey
                          ) {

  # work in data tables
  d <- data.table(data)

  # mother_id
  setnames(d,'mid','mother_id')
  d[, mother_id := paste0(nid,'-',trim(mother_id))]
  d[, child_id  := paste0(mother_id,'-',childs_line)]

  # get P1, P2, P3, average parities for women in the bottom 3 maternal age bins (15-19, 20-24, 25-29), for each NID
  parity <- unique(d[,c('nid','mother_id','mothers_age_group','ceb'),with=FALSE])
  parity <- subset(parity, mothers_age_group %in% c(1,2,3))
  parity <- parity[,.(P=sum(ceb)/.N), by=.(nid,mothers_age_group)]
  parity <- reshape(parity, idvar = "nid", timevar = "mothers_age_group", direction = "wide")
  setnames(parity,old=c('P.1','P.2','P.3'),new=c('P1','P2','P3'))
  d <- merge(d,parity,by='nid',all=TRUE)

  # keeping only data from mothers (CEB>0)
  d <- subset(d, ceb > 0)

  # get percent of ceb by that childs b-day
  d <- d[order(mother_id,child_dob_cmc)]
  d[, birthorder := 1:.N,by=mother_id]
  d[, pct_ceb := birthorder/ceb]
  d[, firstborn := birthorder==1]

  # Cleaning: remove information at mother level if key variables are missing for any children
  d[, na.drop := is.na(birthtointerview_cmc)]
  d[, na.drop := is.na(child_alive)+na.drop]
  d[, na.drop := is.na(interview_date_cmc)+na.drop]
  d[, na.drop := is.na(child_age_at_death_months)+na.drop]
  d$na.drop[d$na.drop>1] <- 1
  d[,na.drop:=max(na.drop),by=mother_id]
  message(sprintf('Dropping %i of %i children (from %i mothers), due to missingness.',
      sum(d$na.drop),nrow(d),length(unique(d$mother_id[d$na.drop==1]))))
  d <- subset(d,na.drop == 0)

  # rename some things
  setnames(d,old='child_age_at_death_months',new='aod')

  # Sort out some year timing
  d[,yrborn    := cmc_yr_convert(interview_date_cmc-birthtointerview_cmc,zero_year=zero_year)]
  d[,yrint     := cmc_yr_convert(interview_date_cmc,zero_year=zero_year)]
  d[,monthint  := round(yrint  * 12,0)]
  d[,monthborn := round(yrborn * 12,0)]

  # calculate mothers age at birth. Round to nearest year.
  d[,mothage_atbirth := round(mothers_age - (yrint - yrborn),0)]
  setnames(d,'mothers_age','mothage_atint')

  # drop if MAB is below 12
  d[,na.drop:=ifelse(mothage_atbirth<10,1,0)]
  d[,na.drop:=max(na.drop),by=mother_id]
  message(sprintf('Dropping %i of %i children (from %i mothers), due to MAB < 10.',
      sum(d$na.drop),nrow(d),length(unique(d$mother_id[d$na.drop==1]))))
  d <- subset(d,na.drop == 0)

  # we are concerned with mortality risk under age 5 for now. Make alive indicator
  # throw a check to make sure the 'alive' indicator was coded as 6000
  stopifnot(get_mode(d$aod)==alive_code)

  # alive, will get used as a censoring indicator
  d[,alive := aod == alive_code]   # if alive at time of survey aod is coded as 6000 as J/A perpped it
  d[,alive := aod >= censor_age]   # Consider all those outside the study period as alive
  d[,died  := !alive]

  # make a child age at time of survey (in months) variable, if died, use AOD
  # This is the random variable T, which may or may not be censored
  d[,childage := monthint-monthborn]      # if alive
  d$childage[!d$alive] <- d$aod[!d$alive] # if died
  d[,childage := childage + 0.0001]

  # how was child age at death recorded (day, month, year)
  d$aod_record <- suppressWarnings(
                    c('d','m','y')[as.numeric(substr(d$child_age_at_death_raw,1,1))] )

  # CD/CEB
  d[, cdceb    := ced/ceb]
  d[, cdceb100 := cdceb*100]

  # return
  return(d)
}







# function to pull all SBH-relevant information for a given mother (mother_id) in the dataset
get_sbh_mother_info <- function(data,          # dataframe that the ID is in, to subet
                                ids = NULL,    # Mother IDs to keep (char vector)
                                form = NULL,          # model formula for prediction
                                variables = NULL # supplied if form = NULL
                                   ){


  # get vars from formula and others we need
  if(!is.null(form) & is.null(variables)){
    form <- strsplit(paste(form)[3],' \\+ | \\* | \\: ')[[1]]

    # remove factor from any variables
    form <- gsub('\\)','',gsub('factor\\(','',form))
    vars <- form[!form %in% c('-1','ab','mothage_atbirth','yrborn','sdi')]

  } else if (is.null(form) & !is.null(variables)){
    vars <- variables[!variables %in% c('mothage_atbirth','yrborn','sdi')] # these get imputed later
  } else {
    stop('variables and form cannot both be NULL')
  }

  # all vars we may need, add them on
  vars <- unique(c(vars,'mothage_atint','yrint','nid','country','ceb','weight'))

  # all mother if null
  if(is.null(ids)) ids <- unique(data$mother_id)

  # keep mother level info only
  data <- data[,c(lapply(.SD,first)), by = mother_id, .SDcols = vars][mother_id %in% ids]

  return(data)

}








# A function to make prediction frame from prepped DHS dataset given mother data from get_sbh_mother_info()
# it will be long by year or month of birth and discrete child age group
expand_prediction_frame <- function(moth_data,               # output data.table from get_sbh_mother_info()
                                    resolution  = 'yearly',  # yearly or monthly output resolution
                                    min_mothage = 10,        # min age a mother can be in prediction frame
                                  #  min_year    = 10,        # min year of birth in prediction frame (20=1990)
                                    abt         = ab_times, # age bins for prediction frame
                                    mothageatintvar = 'mothage_atint'
                                         ){

  # set frame resolution
  stopifnot(resolution %in% c('yearly','monthly'))
  reso <- ifelse(resolution=='yearly',1,1/12)

  # expand mothers ages and dates going back to age <min_mothage>
  moth_data$mato <- moth_data[[mothageatintvar]]
  moth_data <- moth_data[, list(mothage_atbirth = seq(min_mothage,max(mato),reso)), by = names(moth_data)]
  moth_data[,yrborn := yrint - (mato - mothage_atbirth)]
  if(resolution=='yearly') moth_data[,yrborn := round(yrborn,0)]

  # thin out as required
  #moth_data <- subset(moth_data, yrborn >= min_year)

  # make an id for each hypothetical child
  moth_data[, hyp_child_id := paste0(mother_id,'-',1:.N), by=mother_id]

  # expand once more on age bins
  moth_data <- moth_data[, list(ab = abt$ab), by = names(moth_data)]
  moth_data <- merge(moth_data,abt,by='ab')

  # order columns
  cols <- c('mother_id','hyp_child_id','ab')
  setcolorder(moth_data, c(cols, setdiff(names(moth_data), cols)))

  # include some transformations
  dl[, country_ab := as.factor(paste0(country,'_',ab))]

  # return
  return(moth_data)
}








# prediction for all hypothetical children-agebins
predict_Dsurv <- function(new_data,
                          model      = mod,        # glm object that can take predict
                          ci         = TRUE,        # CI on estimates
                          form.fixed = NULL,
                          ndraws     = 250,
                          summarize  = TRUE,
                          ncores     = 1,
                          yrbornvar = 'yrborn'   ) {

  # hazard conditional probability (qx)
  if(class(model)[1] != 'inla'){
    message('Predicting for non-INLA model.')

    # get ci if requested
    if(ci==TRUE) {
      if(any(class(mod) == 'gam')){
        pred <- predict_GAM(m=mod,newdata=new_data,transformation='plogis',draws=ndraws,cores=ncores,summarize=summarize)
        if(summarize == TRUE){
          new_data$haz  <- pred$pred
          new_data$haz_lower   <- pred$lower
          new_data$haz_upper   <- pred$upper
        } else {
          new_data$haz <- pred[,apply(.SD,1,quantile,probs=c(0.5)),.SDcols=colnames(pred)]
        }

      } else { # GLM
        pred <- predict(model,new_data,se.fit=TRUE)$se
        new_data$haz  <- pred$fit
        new_data$se   <- pred$se
        new_data[,haz_lower := plogis( (haz) - 1.96*se)]
        new_data[,haz_upper := plogis( (haz) + 1.96*se)]
        new_data[, haz := plogis(haz)]
      }

    } else {
      new_data$haz <- plogis(predict(model,new_data,se.fit=FALSE))
    }
  } else {
    message('Predicting for INLA model.')
    ff <- as.formula(paste0('~',as.character(form.fixed)[3])) # remove the outcome
    pred <- predictINLA(res=model,newdat=new_data,fixed.formula=ff,summary=TRUE,draws=ndraws,transform='plogis')
    new_data[,haz := pred[,1]]
    new_data[,haz_lower := pred[,2]]
    new_data[,haz_upper := pred[,3]]
  }

  # survival is cumulative hazard in discrete bins
  new_data[,surv := cumprod(1-haz), by=.(mother_id, hyp_child_id)]

  # time of risk associated with each hyp child and age bin
  new_data[,year_entering := floor(get(yrbornvar) + tstart/12) ]

  # subset out any predictions in year of survey as these are hard to compare EEB due to mid-year censoring in the data
  # all predictions end in period before survey year
  new_data[,fyear := floor(yrint),by=nid]
  drop <- new_data$year_entering < new_data$fyear

  # return
  if(summarize) {
    return(new_data[drop,])
  } else {
    return(list(new_data = new_data[drop,], draws = pred[drop,]))
  }
}





## DEPRECIATED, MOVED TO MGCV
# predict out INLA fit, including with random effects (intercept only for now) at new levels
predictINLA <- function(res,draws,fixed.formula,newdat,seed=123445,random.slope=NULL,summary=FALSE,transform=NULL){

  # make samples
  message('Getting posterior samples from INLA ... \n')
  s        <- inla.posterior.sample(result=res,n=draws,seed=seed)

  # get idx and names for fixed and random effects
  message('Extracting Covariate and RE names ... \n')
  fenames  <- res$names.fixed
  revars   <- names(res$summary.random)
  #feidx    <- which(grepl(paste0(fenames,collapse='|'),rownames(s[[1]]$latent)))
  feidx    <- which(!grepl(paste0(c(revars,'Predictor'),collapse='|'),rownames(s[[1]]$latent)))
  stopifnot(length(feidx)==length(fenames))
  stopifnot(all(fenames == rownames(s[[1]]$latent)[feidx]))
  relvls <- reidx <- redraws <- list()
  for(r in revars)  {
    relvls[[r]]   <- res$summary.random[[r]]$ID
    reidx[[r]]    <- which(grepl(paste0(r,'\\:',collapse='|'),rownames(s[[1]]$latent)))
    redraws[[r]]  <- do.call(rbind, lapply(s,function(x){x$latent[reidx[[r]]]}))
    colnames(redraws[[r]]) <- relvls[[r]]
    redraws[[r]] <- data.table(redraws[[r]])
  }

  # pull out the draws for coefficients
  fedraws  <- do.call(rbind, lapply(s,function(x){x$latent[feidx]}))

  # get RE hyperparameter info for levels unknown
  revariances <- (res$summary.hyperpar[paste0('Precision for ',revars),1])^(-1) # TODO later get hyper param draws as well
  for(r in 1:length(revars))  message(paste0('SD for RE <',revars[r],'> is ',round(sqrt(revariances[r]),5)))

  # for each new level in the new data simulate draws
  message('Simulating RE draws for detected new levels ... \n')
  for(r in revars){
    nl <- unique(newdat[[r]])[!unique(newdat[[r]]) %in% relvls[[r]]] # find new levels
    message(paste0(length(nl),' new levels of ',r,' found.'))
    for(lev in nl)
      redraws[[r]][[as.character(lev)]] <- rnorm(draws,0,sqrt(revariances[which(r %in% revars)]))
  }

  # idenify column numbers of random effect draws
  for(r in revars)
    newdat[[paste0(r,'_idx')]] = match(newdat[[r]],colnames(redraws[[r]]))

  # for each row, calculate the fixed effects
  # TODO, generalize fixed eff formula
  message('Making final predictions ... ')
  fmm       <- model.matrix(fixed.formula,newdat)
  fe.dr     <- fmm %*% t(fedraws)

  # for each row, grab the RE if they are in sample, if not, simulate them
  if(length(revars)>0) rmm   <- newdat[,paste0(revars,'_idx'),with=FALSE]

  # random slope
  if(!is.null(random.slope)){
    message('Predicting for random slopes ... ')
    rs.dr <- fmm[,random.slope['fixed']] %*% t(fedraws[,which(fenames==random.slope['fixed'])])
    rs.dr <- rs.dr * t(redraws[[random.slope['random']]][,rmm[[paste0(random.slope['random'],'_idx')]],with=FALSE]) # !! SOME ERROR HERE CHECK IT
    revars <- revars[! revars %in% random.slope['random']]
  } else {
    rs.dr <- matrix(0,nrow=nrow(newdat),ncol=draws)
  }

  # simulate random intercepts
  re.dr <- matrix(0,nrow=nrow(newdat),ncol=draws)
  for(r in revars) re.dr <- re.dr + t(redraws[[r]][,rmm[[paste0(r,'_idx')]],with=FALSE])

  # combine fixed and random effects
  pred.dr <- fe.dr + re.dr + rs.dr

  # tranform in draw space if requested
  if(!is.null(transform)){
     message(paste0('Transforming using ',transform))
     pred.dr <- eval(parse(text=paste0(transform,'(pred.dr)')))
  }

  # data.table
  pred.dr <- data.table(pred.dr)

  if(summary){
    message('Summarizing draws ... \n')
    pred.dr <-t(pred.dr[,apply(.SD,1,quantile,probs=c(0.5,0.025,0.975)),.SDcols=colnames(pred.dr)])
  }

  return(pred.dr)
}







#### function to predict for new levels
# hold out any data with new levels
predict_GAM <- function(m ,                           # gam model object with s(...,bs='re')
                        newdata,                      # newdata frame
                        draws = 250,                 # prediction draws
                        transformation = NULL,         # link
                        parallel = TRUE,
                        summarize = TRUE,
                        cores  = 1
                        ){

  message('Predicting for mgcv::gam or mgcv::bam object with potential s(..,bs="re") style random effects ... \n ')
  message('WARNING: Make sure random effect variable names are unique or else there may be grep conflicts ... \n ')

  # requred library function for drawing from a multivariate normal distribution
  require(MASS)
  require(data.table)
  require(parallel)

  # get cores to use
  message(sprintf('Using %i cores ... \n ',cores))

  nd <- data.table(copy(newdata))
  nd$idx <- 1:nrow(nd)

  # identify the random effect(s) from the model formula
  f <- m$formula
  fparts  <- strsplit(as.character(f)[3], ' \\+ ')[[1]]
  random  <- fparts[grep('bs = \"re\"',fparts)]
  message(paste0(length(random),' random effects detected: '))
  if(length(random)!=0){
    for(i in 1:length(random))
      random[i] <- gsub(', bs = \"re\"\\)','',gsub('s\\(','',random[[i]]))
      message(paste(random,collapse = ', '))

    # split newdata into new levels and old levels
    for(r in random) nd[[paste0(r,'_newlevel')]] <- ! as.character(newdata[[r]]) %in%  as.character(unique(m$model[,r]))
    nd <- nd[,oos_levels := rowSums(.SD) ,.SDcols=paste0(random,'_newlevel')]
    nd1 <- nd[oos_levels == 0]
    nd2 <- nd2_orig <- nd[oos_levels >  0]
    message(sprintf('Detected new levels in %i of the random effects ... \n ',max(nd$oos_levels)))
  } else {
    message('No random effects detected ... \n\n ')
    nd1 <- nd
    nd2 <- nd[0]
  }

  #####
  # in sample all levels is easy, so do that first
  if(nrow(nd1) > 0 ){
    message('Predicting for rows with in-sample random effects levels ... \n ')

    message(' ... predicting out lp matrix ')
    cl     <- makeCluster(cores)
    X <- predict(m, nd1, type="lpmatrix") # get bases
    stopCluster(cl)

    message(' ... drawing betas ')
    b1_d <- mvrnorm(n = draws, mu = m$coefficients , Sigma = m$Vp) # draw betas
    if(cores > 1){ # in parallel helps
      message(' ... IS matrix multiplication for rows with IS levels in parallel ')
      cl     <- makeCluster(cores)
      idx    <- splitIndices(nrow(X), length(cl))
      Xlist  <- lapply(idx, function(r) X[r,,drop=FALSE])
      p1_d   <- do.call(rbind, clusterApply(cl, Xlist, get('%*%'), t(b1_d)))
      stopCluster(cl)
    } else {
      message(' ... IS matrix multiplication for rows with IS levels NOT in parallel ')
      p1_d <- X %*% t(b1_d)
    }
    message(' ... apply transformation at the draw level ')
    if(!is.null(transformation)) p1_d <- match.fun(transformation)(p1_d) # transform drawl level preds
    p1_d <- data.table(p1_d) # summarize
    if(summarize) {
      message(' ... summarizing draws into 0.025, 0.50, and 0.975 quantiles ')
      p1   <- data.table(t(p1_d[,apply(.SD,1,quantile,probs=c(0.5,0.025,0.975)),.SDcols=colnames(p1_d)]))
      colnames(p1) <- c('pred','lower','upper')
    } else {
      p1 <- p1_d
    }
    p1$idx <- nd1$idx
  } else {
    p1 <- data.table() # if all nd rows have oos levels
  }
  #####
  # loop through random variables and for oos levels draw from their estimated distribution
  if(nrow(nd2) > 0){
    message('Predicting for rows with new levels of random effects by drawing from estimated RE distribution ... \n ')

    message(' ... make fake levels to trick predict.gam into giving us lpmatrix ')

    for(r in random){ # TODO CHECK ON THIS FACTOR REPLACING IS WEIRD SO MAKE SURE ITS DOING WHAT YOU THINK
      if(all(nd2[[paste0(r,'_newlevel')]] == TRUE)){
        nd2[[r]] = as.factor(data.table(m$model)[[r]][1])
      } else {
        nd2[[r]][nd2[[paste0(r,'_newlevel')]] == TRUE] = as.factor(data.table(m$model)[[r]][1])
      }
    }

    # split out in sample and out of sample RE bases
    message(' ... predicting out lp matrix ')
    cl     <- makeCluster(cores)
    X <- predict(m,nd2,type="lpmatrix",cluster=cl)
    stopCluster(cl)

    # turn X indicators into all zeros for each missing RE
    message(' ... X to zero ')
    for(r in random)
      for(v in grep(paste0(r,collapse='|'),colnames(X)))
        X[nd2[[paste0(r,'_newlevel')]] == TRUE,v]=0

    # IS predictions at draw level
    message(' ... drawing betas ')
    b1_d <- mvrnorm(n = draws, mu = m$coefficients , Sigma = m$Vp) # draw betas
    if(cores > 1){ # in parallel helps
      message(' ... IS matrix multiplication for rows with OOS levels in parallel ')
      cl     <- makeCluster(cores)
      idx    <- splitIndices(nrow(X), length(cl))
      Xlist  <- lapply(idx, function(r) X[r,,drop=FALSE])
      p2_d   <- do.call(rbind, clusterApply(cl, Xlist, get('%*%'), t(b1_d)))
      stopCluster(cl)
    } else {
      message(' ... IS matrix multiplication for rows with OOS levels NOT in parallel ')
      p2_d <- X %*% t(b1_d)
    }

    # for any missing RE now draw those from their estimated distribution and add them on to the predictive draws
    message(' ... grab estimated variances from random effects and use them to draw rnorm ')
    for(r in random) {
      gvcmp <- mgcv::gam.vcomp(m)
      re_sd <- gvcmp[grepl(r,rownames(gvcmp)),1] # TODO, use lower and upper from this as well
      rows  <- which(nd2[[paste0(r,'_newlevel')]] == TRUE)
      p2_d[rows,]  <- p2_d[rows,] + rnorm(length(rows)*draws,0,re_sd)
    }

    # apply transform
    message(' ... apply transformation at the draw level ')
    if(!is.null(transformation)) p2_d <- match.fun(transformation)(p2_d) # transform drawl level preds
    p2_d <- data.table(p2_d) # summarize
    if(summarize) {
      message(' ... summarizing draws to 0.025%, 0.50%, and 0.975% ')
      p2 <- data.table(t(p2_d[,apply(.SD,1,quantile,probs=c(0.5,0.025,0.975)),.SDcols=colnames(p2_d)]))
      colnames(p2) <- c('pred','lower','upper')
    } else {
      p2 <- p2_d
    }
    p2$idx <- nd2$idx

    # combine p1 and p2
    message(' ... combine IS and OOS data chunks ')
    cpred <- rbind(p1,p2)
  } else {
    message('No out of sample levels in Random effects were detected ... \n ')
    cpred <- p1
  }

  #  make sure we return it in the same oder we got it from newdata
  message('Order data ... \n')
  setorder(cpred, idx)
  cpred$idx <- NULL

  return(cpred)
}







# get IHME LOC INFO
merge_ihme_loc <- function(d      # d must have a variable country that is ISO3
                              ){
  d <- merge(d,fread("/share/geospatial/royburst/sbh_cbh/data/ihme_lc_info.csv"),
              by='country',all.x=TRUE)
  return(d)
}







# load MAP distributions
load_MAP_distributions <- function(gbdregions){

  load(sprintf('%s/temp/geospatial/U5M_global/sbh_cbh/data/MAP.rdata',j()))
  load('/share/geospatial/royburst/sbh_cbh/data/MAP.rdata')

  map <- subset(data.table(distributions$ceb),gbdregion %in% gbdregions)
  map[,agegroup:=as.numeric(substr(agegroup, 1, 2))]

  return(map)
}





## Make time stamp in standardized format.
make_time_stamp <- function(time_stamp=TRUE) {

    run_date <- gsub("-","_",Sys.time())
    run_date <- gsub(":","_",run_date)
    run_date <- gsub(" ","_",run_date)

    if(time_stamp==FALSE) run_date <- 'scratch'

    return(run_date)

  }




## Load config data into environment
load_config <- function(dir){
  config <- fread(sprintf('%s/config.csv',dir))[,2:3]
  message(' ... LOADING CONFIG VARIABLES INTO NAMESPACE:')
  for(param in config[, V1]) {
    tmp <- config[V1==param, V2]
    tmptest <- suppressWarnings(try(as.numeric(tmp)))
    tmp <- ifelse(is.na(tmptest),tmp,tmptest)
    message(sprintf(' ... ... %s ::: %s',param,tmp))
    assign(param, tmp, envir=globalenv())
  }
}





# make a country dir in the standard format
# ./output/<run_date>/<country>/. Country is the one being held out
make_output_dirs <- function(root     = '/share/geospatial/royburst/sbh_cbh/output',
                             run_date,
                             country  = NULL) {

  # make the run_date directory
  suppressWarnings(
    dir.create(sprintf('%s/%s', root, run_date))
  )

  # make a country directory within the run_date directory
  if(!is.null(country))
    suppressWarnings(
      dir.create(sprintf('%s/%s/%s', root, run_date, country))
    )
}





# center scale
centscale <- function(var,varname,centscale=cs,reverse=FALSE){
  if(varname %in% centscale$name){
    if(!reverse)
      res <- round((var - centscale$mean[centscale$name==varname])/ centscale$sd[centscale$name==varname] ,2)
    else
      res <-round((var * centscale$sd[centscale$name==varname] ) + centscale$mean[centscale$name==varname],2)
  } else {
    message(sprintf('Varname %s not found in cs table, returning same variabe',varname))
    res <- var
  }
  return(res)
}



# get a table for center scaling
get_cs_table <- function(data,variables,skipvars){
  cs <- data.table(name=character(),mean=numeric(),sd=numeric())
  for(v in variables){
    if(class(data[[v]]) == 'numeric' | class(data[[v]]) == 'integer' & ! v %in% skipvars){
      message(sprintf('Center Scaling:      %s ', v))
      data[[paste0(v,'_orig')]]=data[[v]]
      sc     <- scale(data[[v]],center=TRUE,scale=TRUE)
      data[[v]] <- sc
      cs <- rbind(cs,data.table(name=v,mean=attr(sc,"scaled:center"),sd=attr(sc,"scaled:scale")))
    } else if(v %in% skipvars){
      message(sprintf('Skipping Variable:      %s ', v))
      cs <- rbind(cs,data.table(name=v,mean=0,sd=1))
    } else {
      stop(sprintf('Found FE variable %s which is not numeric or in skipvars list..... ', v))
    }
  }
  return(cs)
}





# parse out all the fixed effects variables from the gam formula
get_fe_from_gamformula <- function(gamformula){

  # first split out components
  form <- gamformula
  form <- gsub('\"','',form)
  form <- gsub("\'",'',form)
  form <- unlist(strsplit(form, ' \\+ '))
  form <- unlist(strsplit(form, '\\*'))

  # get random effects, knock them out, since we only want fe right now
  re <- form[grepl('\\=re|\\= re',form)]
  message(sprintf('Random effects found in formula: %s ...\n',paste(re,collapse='  &  ')))
  form <- form[!grepl('\\=re|\\= re',form)] # remove REs

  # remove nusances
  form <- form[!grepl('\\~',form)]  # remove the equals bit
  #form <- form[grepl('s\\(',form)] # only take smooths

  # clean up so we only get a list of FEs
  form <- unlist(strsplit(gsub('s\\(','',form),', '))
  form <- unlist(strsplit(gsub('te\\(','',form),', '))
  form <- form[!grepl('k \\=|k\\=',form)]
  vars <- form[!grepl('by = factor\\(ab\\)|by=factor\\(ab\\)|by =factor\\(ab\\)|by= factor\\(ab\\)',form)]
  vars <- form[!grepl('factor\\(ab\\)',form)]

  message(sprintf('Fixed effects found in formula: %s ...\n',paste(vars,collapse='  &  ')))
  return(vars)

}






# make aggregated estimates across period
aggregate_period_estimates <- function(est,   # predictions at the hyp child level
                                       map,   # MAP distributions
                                       use_svy_weight  = FALSE, # look for a variable called 'svywgt' and use it for aggregation
                                       ci = TRUE,
                                       TESTNOMAP = FALSE,
                                       yrbornvar = 'yrborn',
                                       cebvar    = 'ceb',
                                       mothageatintvar = 'mothage_atint',
                                       yrs_to_agg = 1,
                                       return_individual = FALSE) {

  if(class(est) == 'list') { # this means it was unsummarized, with a new_data object and a draws object
    message('Detected draws level data, so aggregating by draws .. ')
    draws <- est$draws
    est   <- est$new_data
    draws$idx <- 1:nrow(draws)
    est$idx   <- 1:nrow(est)
    drawlevel <- TRUE
  } else {
    message('Detected summary level data, so aggregating by summary quantiles .. ')
    drawlevel <- FALSE
  }

  # turn mothers age at survey into two year interval to match with MAP distribution
  est$mato <- est[[mothageatintvar]]
  est$mato[est$mato%%2==1] <- est$mato[est$mato%%2==1]-1
  est$mato[est$mato<=16]   <- 15

  # merge
  est <- merge(est,map,by.x=c('mato',cebvar),by.y=c('agegroup','ceb'),all.x=TRUE)

  # get weight t.X based on time birth is prior to survey
  est[, yrsprior := -1*(round(get(yrbornvar)-yrint,0))-1 ] # TDOfix this to be more exact (issues from rounding errors in yrborn earlier.)
  #est[, yrsprior2 := -1*(floor(get(yrbornvar)-yrint)) ] # fix this to be more exact
  #est[, yrsprior3 := floor(-1*(get(yrbornvar)-yrint)) ] # fix this to be more exact

  est <- subset(est, yrsprior <= 24) # This can cause comparison problems: make sure cbh is subset to 24 yrs prior as well

  # get the expected number entering each bin
  est$wgt <- NA
  for(i in 0:24) {
    if(TESTNOMAP == FALSE) {
      est$wgt[est$yrsprior==i] <- est[[paste0('t.',i)]][est$yrsprior==i]
    } else  { # if testing no MAP dist, set all weights to one
      est$wgt[est$yrsprior==i] <- 1
    }
    est[[paste0('t.',i)]]    <- NULL
  }
  if(TESTNOMAP == FALSE) {
    est$EEB <- est$wgt*est[[cebvar]]*((est$surv)/(1-est$haz))
  } else {
    est$EEB <- est$wgt*((est$surv)/(1-est$haz))
  }

  if(return_individual) {
    return(est)
  } else { # aggregate

    # make year bins
    yrmap <- data.table(year_entering = seq(from = min(est$year_entering), to = max(est$year_entering)),
                        yrmap=rep(seq(from = min(est$year_entering), to = max(est$year_entering), by = yrs_to_agg),each=yrs_to_agg) + .5*yrs_to_agg)
    yrmap <- subset(yrmap, abs(year_entering - yrmap)<yrs_to_agg )
    est   <- merge(est, yrmap, by='year_entering')

    # WHy is pnn2 double counted?
    # head(est[,.(eeb=.N),by=.(ab,year_entering)],30) # 2 are entering each year is why!

    # make new weights based on EEB and svywgt, if its available
    if(use_svy_weight == FALSE) est[, svywgt := 1]
    est[, combined_weight := EEB*svywgt]


    # aggregate over age bin and year entering.
    if(!ci){
      est <- est[,.(haz = weighted.mean(haz,wgt),EEB=sum(EEB)), by = .(ab,yrmap)]
    } else {
      if(drawlevel){ # DO BY DRAW
        message('Aggregating by draws')
        drawcols  <- colnames(draws)[!colnames(draws) %in% 'idx']
        draws     <- merge(draws, est[,c('EEB','ab','yrmap','idx','combined_weight'),with=FALSE], by = 'idx', all.y = TRUE)
        draws_agg <- draws[, lapply(.SD, weighted.mean, combined_weight), .SDcols = drawcols, by = .(ab,yrmap)]
        eeb_agg   <- draws[, .(EEB = sum(EEB)), by = .(ab,yrmap)]
        summ      <- data.table(t(draws_agg[,apply(.SD,1,quantile,probs=c(0.5,0.025,0.975)),.SDcols=drawcols]))
        colnames(summ) <- c('haz','haz_lower','haz_upper')
        res       <- cbind(eeb_agg,summ)
        setnames(draws_agg, 'yrmap', 'year_entering')
      } else { # THIS IS A SHORTCUT AND SHOULD NOT GENERALLY BE DONE
        message('Aggregating by quantiles ... do not do for final runs!')
        res <- est[,.(haz       = weighted.mean(haz,EEB),
                      haz_lower = weighted.mean(haz_lower,EEB),
                      haz_upper = weighted.mean(haz_upper,EEB),
                      EEB=sum(EEB)), by = .(ab,yrmap)]
      }
    }

    # return
    setnames(res, 'yrmap', 'year_entering')

    if(drawlevel){
      return(list(res=res,draws=draws_agg))
    } else {
      return(res)
    }
  }
}






# get pct ceb from MAP for predictions
get_pct_ceb <- function(newdata = nd, mapp = map, mothageatintvar = 'mothage_atint_orig', cebvar='ceb_orig' ,    yrbornvar = 'yrborn_orig' ){

  # turn mothers age at survey into two year interval to match with MAP distribution
  newdata$mato <- newdata[[mothageatintvar]]
  newdata$mato[newdata$mato%%2==1] <- newdata$mato[newdata$mato%%2==1]-1
  newdata$mato[newdata$mato<=16]   <- 15

  # merge
  newdata <- merge(newdata,mapp,by.x=c('mato',cebvar),by.y=c('agegroup','ceb'),all.x=TRUE)
  newdata[, yrsprior := (round(yrint-get(yrbornvar),0)) ]

  # make cumulative percent born
  for(i in 23:0) newdata[[paste0('t.',i)]]<- newdata[[paste0('t.',i)]]+newdata[[paste0('t.',i+1)]]
  # for older women, add on any missing part of the distribution
  for(i in 24:0) newdata[[paste0('t.',i)]] <-  newdata[[paste0('t.',i)]] + (1-newdata[[paste0('t.0')]])
  newdata$yrsprior[newdata$yrsprior>24] <- 24
  newdata$pct_ceb_orig=NA
  for(i in 24:0) {
    newdata$pct_ceb_orig[newdata$yrsprior==i] <- newdata[[paste0('t.',i)]][newdata$yrsprior==i]
    newdata[[paste0('t.',i)]] <- NULL
  }
  newdata$ceb.group <- newdata$mato <- newdata$gbdregion <- NULL

  return(newdata)
}





# get raw hazards from test and train data
# collapse on age bin and year entering the bin to get an empirical period estimate for that bin
tabulate_raw <- function(earliest_year, dl_train, dl_test, by = base_year, yrs_to_agg = 1, oos, use_svy_weight = FALSE){

  for(tt in c('test','train')){
    message(paste0(' ... ... ... ',tt))
    aggdat <- subset(get(paste0('dl_',tt)), year_entering < fyrint)

    # make year bins
    yrmap <- data.table(year_entering = seq(from = min(aggdat$year_entering), to = max(aggdat$year_entering)),
                        yrmap=rep(seq(from = min(aggdat$year_entering), to = max(aggdat$year_entering), by = yrs_to_agg),each=yrs_to_agg) + .5*yrs_to_agg)
    yrmap <- subset(yrmap, abs(year_entering - yrmap)<yrs_to_agg )
    aggdat <- merge(aggdat, yrmap, by='year_entering')

    if(tt == 'test'){ # this only works for one survey, because qweights are diffferent between surveys
      if(use_svy_weight == FALSE) aggdat[, weight := 1]
      aggdat[, weight   := weight/sum(weight,na.rm=TRUE)]
    } else{
      aggdat[, weight := 1]
    }

    if(oos) aggdat <-  na.omit(aggdat[,.(haz=weighted.mean(died, weight),entering=sum(dum)),by=.(ab,yrmap, country)])
    if(!oos) aggdat <- na.omit(aggdat[,.(haz=weighted.mean(died, weight),entering=sum(dum)),by=.(ab,yrmap)])
    aggdat <- subset(aggdat, yrmap >= (earliest_year-base_year))
    aggdat[,train := tt=='train']
    aggdat[, year := yrmap+by]
    setnames(aggdat, 'yrmap', 'year_entering')
    assign(paste0('aggdat_',tt),aggdat)
  }
  return(rbind(aggdat_test,aggdat_train))
}





# get q5 from a dataset with all exclusive ages underneath it
get_q5 <- function(data, hazvars){
  q5 <- copy(data)
  q5[, a := as.numeric(substr(as.character(ab),1,1))]
  q5 <- q5[order(train,year,a)]  # order by age bin and period
  q5 <- q5[,n := sum(.N),by=.(train,year)] # remove any periods without all age bin
  q5 <- subset(q5, n == nrow(ab_times))
  for(a in 2:nrow(ab_times))
    for(hazvar in hazvars)
      q5[[hazvar]][q5$a==a]   = 1-(1-q5[[hazvar]][q5$a==a])   * (1-q5[[hazvar]][q5$a==(a-1)])

  q5 <- subset(q5, a==n) # keep only last row
  q5$a <- q5$n <- NULL
  q5$ab <- '8_5q0'
  return(q5)
}





# calculate some oos pv metrics, using the comp dataset (comparison of est and dat)
get_oos_pv_summary_metrics <- function(data){
  pv <- subset(data, train==FALSE)
  pv <- pv[, q_error := haz_est - haz_dat]
  pv <- pv[, N_error := EEB - entering]
  pv <- pv[, .( q_mean_est = mean(haz_est),
                q_mean_dat = mean(haz_dat),
                q_rmse     = sqrt( mean( q_error^2 ) ),
                q_me       = mean(q_error),
                q_95cov    = mean(haz_dat>haz_lower & haz_dat<haz_upper),
                N_mean_est = mean( EEB ),
                N_mean_dat = mean( entering ),
                N_rmse     = sqrt( mean( N_error^2 ) ),
                N_me       = mean(N_error)),
            by = .(ab)  ]

 return(pv)
}




dashplots <- function(filename,oos=FALSE,form=NULL){
  #ggplot(aggres,aes(x=year,y=haz,colour=ab,group=ab))+

  # model versus the data facetted by bin
  if(! 'country' %in% colnames(aggdat)) aggdat$country = 'XXX'
  g1=ggplot(aggres,aes(x=year,y=haz,colour=ab,group=ab,fill=ab))+
    geom_line(data=aggdat[train==TRUE],color='black',aes(group=country),alpha=.1)+
    geom_line(aes(y=haz_lower),size=0.125,linetype=2,lineend='round')+
    geom_line(aes(y=haz),size=0.5,linetype=1,lineend='round')+
    geom_line(aes(y=haz_upper),size=0.125,linetype=2,lineend='round')+
    ggtitle('95%CI of Hazard (q_x) Estimates over Time')+facet_grid(~ab)+
    scale_x_continuous(breaks = seq(1985,2020,by=10))+
    theme(strip.background =element_rect(fill="white"))+
    theme_bw()
  if(oos) g1 <- g1 + geom_line(data=aggdat[train==FALSE],color='red',alpha=.7)

  if(oos) comp = comp[train==FALSE] else comp = comp[train==TRUE]
  # scatter of estimates
  g2=ggplot(comp,aes(y=haz_est,x=haz_dat,col=ab))+
    geom_point(aes(size=entering))+theme_bw()+
    geom_errorbar(aes(ymin=haz_lower,ymax=haz_upper))+
    geom_abline(intercept=0,slope=1,col='red')+
    xlab('DHS Observed Probability of Death')+
    ylab('Modeled Probability of Death')+
    ggtitle('Predicted vs Observed q_x')+
    theme(legend.position="none")


  # scatter of number entering each bin
  g3=ggplot(comp,aes(y=EEB,x=entering,col=ab))+
    geom_point(size=2)+theme_bw()+
    geom_abline(intercept=0,slope=1,col='red')+
    xlab('Actual Number Entering Bin-Year')+
    ylab('Expected Number Entering From MAP/Model')+
    ggtitle('Predicted vs Observed SS (# Entering Bin)')+
    theme(legend.position="none")

  # expected number of deaths
  g4=ggplot(comp,aes(y=EEB*haz_est,x=entering*haz_dat,col=ab))+
  #  geom_errorbar(aes(ymin=haz_lower*EEB,ymax=haz_upper*EEB))+
    geom_point(size=2)+theme_bw()+
    geom_abline(intercept=0,slope=1,col='red')+
    xlab('Observed Deaths in each year-bin')+
    ylab('Predicted Deaths in each year-bin')+
    ggtitle('Predicted vs Observed Deaths')+
    theme(legend.position="none")

  gs <- list(g1,g2,g3,g4)

  lay <- rbind(c(1,1,1,1,1,1,1),
               c(1,1,1,1,1,1,1),
               c(2,2,2,3,3,3,3),
               c(2,2,2,3,3,3,3),
               c(2,2,2,4,4,4,4),
               c(2,2,2,4,4,4,4))

  if(length(unique(dl$nid))==1){
    title <- sprintf('NID: %s \n ISO: %s | Svy: %s | Year: %s | OOS: %s',
                    paste0(unique(dl$nid),collapse=', '),dl$country[1],dl$survey[1],dl$year[1],oos)
  } else {
    if(oos == TRUE)
    title <- sprintf('Test NIDs: %s \n Train NIDs: %s \n Countries in Training: %s \n Country in Test: %s',
                    paste0(unique(dl_test$nid),collapse=', '),paste0(unique(dl_train$nid),collapse=', '),paste(unique(dl_train$country),collapse=', '), unique(dl_test$country))
    if(oos == FALSE)
    title <- sprintf(' NIDs: %s \n Country: %s',paste0(unique(dl$nid),collapse=', '),paste(unique(dl$country),collapse=', '))
  }
  if(!is.null(form)) title <-sprintf('%s\n Formula: %s',title,as.character(form)[3])

  g <- grid.arrange(grobs = gs, layout_matrix = lay,
                    top   = textGrob(title,gp=gpar(fontsize=12,font=3)))


  if(!is.null(filename)){
    pdf(filename,width=20,height=15)
    plot(g)
    dev.off()
  } else {
    return(g)
  }

}







q5plot <- function(filename,q5dat,brassdat=NULL,ihmedat=NULL,c=cntry){
  g1=ggplot(q5dat,aes(x=year,y=haz_est))+
  #  geom_line(data=aggdat[train==TRUE],color='black',aes(group=country),alpha=.1)+
    geom_line(aes(y=haz_lower),size=0.125,linetype=2,lineend='round',color='black')+
    geom_line(aes(y=haz_est),  size=1,    linetype=1,lineend='round',color='black')+
    geom_line(aes(y=haz_upper),size=0.125,linetype=2,lineend='round',color='black')+
    geom_line(data=q5dat[train==FALSE],aes(y=haz_dat),color='red',alpha=.7) +
    ggtitle(sprintf('5q0 Estimates Time trend (with 95pct CI) \n red=truth, black=model, blue=brass, green=IHME \n %s',c)) +
    theme_bw()+
    scale_colour_manual(name = '', guide = 'legend',
         values =c('#184231'='#184231','red'='red'), labels = c('OOS estimate','raw (hidden) data'))

  if(!is.null(brassdat)){
    g1 = g1 + geom_line(data=brassdat,aes(y=q5,x=time),color='blue',alpha=.7) +
         geom_point(data=brassdat[3,],aes(y=q5,x=time),color='blue',alpha=.7) # for 25-29 yo
  }
  if(!is.null(ihmedat)){
    g1 = g1 + geom_line(data=ihmedat,aes(y=q5,x=time),color='green',alpha=.7)
  }

  pdf(filename,width=8,height=8)
  plot(g1)
  dev.off()
}


# make qsub to launch a job on IHME's sge
make_qsub     <- function(rtdir,
                          rd,
                          cores,
                          memory,
                          proj,
                          coderepo,
                          codescript,
                          geo_nodes,
                          cntry,
                          adl_args = '',
                          singularity=FALSE) {



  # sort directories
  logloc <- sprintf('%s/%s/%s',rtdir,rd,cntry)
  dir.create(logloc)

  ## do we want to submit to geo nodes? if so, there are a few tweaks
  if(geo_nodes == TRUE){
    shell     <- "r_shell_geos.sh"   ## for the correct path to R
    proj      <- "proj_geo_nodes"    ## correct proj for geos nodes
    node.flag <- "-l geos_node=TRUE" ## send the job to geo nodes
    shell     <- sprintf('%s/shell_geos.sh',coderepo)
  } else {
    proj      <- "proj_geospatial"
    node.flag <- ""
    shell     <- sprintf('%s/shell_prod.sh',coderepo)
  }
  if(singularity == TRUE) shell <- sprintf('%s/shell_singularity10.sh',coderepo)

  # any other args.
  adl_args = paste(adl_args,collapse=' ')

  # code path to run
  code <- sprintf('%s/%s',coderepo,codescript)

  # define qsub
  qsub <- paste0('qsub -j y -o ',logloc,' -cwd -l mem_free=',memory,'G -pe multi_slot ',cores,
                 ' -P ',proj,' ',node.flag,' -N job_',cntry,' ',shell,' ',code,' ',rd,' ',cntry,' ',adl_args,' fin')

  return(qsub)
}
