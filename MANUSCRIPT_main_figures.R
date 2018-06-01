
## #############################################################
## Roy Burstein
## 2018
## Pull results and make figures for use in the manuscript.
## ##############################################################


## --------------------------------------------------------
## SETUP
## load libraries
library(data.table)
library(ggplot2)
library(Hmisc)
library(boot)
library(gtools)

# paths
codepath <- '/share/geospatial/royburst/code/sbh_cbh'     # where the code at
root     <- '/share/geospatial/royburst/sbh_cbh/output'   # dir for outputs

## Make sure code repo being accessed on the cluster is pulled and up to date
system(sprintf('cd %s\ngit pull origin master',codepath))

## load custom functions
setwd(codepath)
source('utils.R')

## get the run date identifying the model runs to use
rds          <- c('2018_03_08_09_56_20',
                  '2018_05_15_13_24_14',
                  '2018_05_15_13_23_29',
                  '2018_05_15_13_22_59',
                  '2018_05_15_13_22_34')

## labels for each model run date
labs <- c('FULL','INDIV','TREND','ADD','INT')

## Of these model runs, note which one is the final full run, as some are comparison runs
FULLMODELRD  <-  '2018_03_08_09_56_20'

## set some paths which depend on the model run
rd_dirs      <- paste0(root,'/',rds)
finaloutputs <- paste0('/share/geospatial/royburst/sbh_cbh/plots/',FULLMODELRD)
## --------------------------------------------------------


## --------------------------------------------------------
## PULL IN ALL THE DATA NEEDED FROM ALL COUNTRIES
## --------------------------------------------------------
## get list of countries as the directories from this run_date output
countries <- unique(list.dirs(path = rd_dirs, full.names = FALSE, recursive = FALSE))

## for each directory, load info needed to make comparison
comps <- list() # each slot for the comparison file for that run date
for(rd in 1:length(rds)){

  message(rd)
  comp <- data.table()

  for(c in countries){
    message(c)

    # load the comparison data, for ages and for q5
    if(file.exists(sprintf('%s/%s/comp.csv',rd_dirs[rd],c))){
      tmp1   <- fread(sprintf('%s/%s/comp.csv',rd_dirs[rd],c))[train==FALSE]
      tmp2   <- fread(sprintf('%s/%s/q5.csv',rd_dirs[rd],c))[train==FALSE]
      tmp    <- rbind(tmp1,tmp2)

      # load ancillary identifying information
      nids       <- fread(sprintf('%s/%s/nids.csv',rd_dirs[rd],c))[train==FALSE]
      tmp$nid    <- as.numeric(nids$nid)
      tmp$svyyr  <- as.numeric(nids$year)
      tmp$region <- as.character(nids$MAP_gbdregion)

      # keep only estimates that are 15 years since survey (temporal scope of this study)
      tmp <- subset(tmp, svyyr-year < 15)

      # rbind them together
      comp <- rbind(comp,tmp)
    }
  }
  comp[,yrsprior := svyyr - year]
  comps[[rds[rd]]] <- comp
}
## --------------------------------------------------------



## --------------------------------------------------------
## ALSO PULL IN AD1 LEVEL COMPARISONS
## --------------------------------------------------------

comps_ad1 <- list() # each slot for the comp file for that run date

## for each directory, load info needed to make comparison
for(rd in 1:length(rds)){
  compad1 <- data.table()
  for(c in countries){
    # load the comparison data, for ages and for q5
    if(file.exists(sprintf('%s/%s/cbh_sbh_compare_ad1_level.csv',rd_dirs[rd],c))){
      ad1   <- fread(sprintf('%s/%s/cbh_sbh_compare_ad1_level.csv',rd_dirs[rd],c))

      if("GAUL_CODE" %in% colnames(ad1)){
        message(sprintf('%s %s',rds[rd],c))
        # load ancillary identifying information
        nids       <- fread(sprintf('%s/%s/nids.csv',rd_dirs[rd],c))[train==FALSE]
        ad1$svyyr  <- as.numeric(nids$year)
        ad1$region <- as.character(nids$MAP_gbdregion)

        # make this data look like the comps data (wide) so we can later use the same
        ad1   <- subset(ad1, nid == nids$nid ) # keep only test id
        ad1   <- subset(ad1, !is.na(GAUL_CODE) ) # ONLY keep identifiable ad1s
        ad1[, year  := year_entering + 0.5]
        ad1[, train := FALSE]
        ad1[, c('V1','age_bin','year_entering') := NULL] # remove some thigns we won't need

        ad1_c <- subset(ad1, source == 'raw cbh')
        setnames(ad1_c, c('N','mr','died'),c('entering','haz_dat','died_dat'))
        ad1_s <- subset(ad1, source == 'indirect sbh')
        setnames(ad1_s, c('N','mr','died'),c('EEB','haz_est','died_est'))
        ad1 <- merge(ad1_c, ad1_s, by=c('ab','GAUL_CODE','NAME','year','nid','svyyr','region','train'))
        ad1[, c('source.y','source.x') := NULL]

        # keep only estimates that are 15 years since survey
        ad1 <- subset(ad1, svyyr-year < 15)
        # rbind them together
        compad1 <- rbind(compad1,ad1)

      } else {
        # this warning should not show up if everything ran and didnt fail
        message(sprintf('WARNING: %s %s skipped because GAUL_CODE not in dataset. make sure you re-ran it',rds[rd],c))
      }

    }
  }
  compad1[,yrsprior := svyyr - year]
  comps_ad1[[rds[rd]]] <- compad1
}
## --------------------------------------------------------



## --------------------------------------------------------
## Subnational, In order to get stable estimates of Ad1, need to aggregate them into 3 yr bins
for(rd in rds){
  comps_ad1[[rd]][yrsprior>10,yrsprior:=12.5]
  comps_ad1[[rd]][yrsprior>5&yrsprior<10,yrsprior:=7.5]
  comps_ad1[[rd]][yrsprior>0&yrsprior<5,yrsprior:=2.5]

  # aggregate adm1 results to the new time bins
  comps_ad1[[rd]] <- comps_ad1[[rd]][, .(
    entering   = sum(entering),
    EEB        = sum(EEB),
    died_dat   = sum(died_dat),
    died_est   = sum(died_est),
    haz_est    = sum(died_est)/sum(EEB),
    haz_dat    = sum(died_dat)/sum(entering)),
    by = .(ab, GAUL_CODE, NAME, nid, svyyr, region, train, yrsprior)]

  comps_ad1[[rd]][, year := yrsprior]
  comps_ad1[[rd]][, loess_haz_dat := haz_dat]   # no smoothing since we aggregated.

  # q5 for comparison
  q5 <- comps_ad1[[rd]][, .(loess_haz_dat  = 1-prod(1-haz_dat),
                            haz_est        = 1-prod(1-haz_est),
                            haz_dat        = 1-prod(1-haz_dat),
                            entering       = sum(entering),
                            EEB            = sum(EEB),
                            died_dat       = sum(died_dat)),
             by = .(NAME,GAUL_CODE,svyyr,yrsprior,year,nid,region)]
  q5[,ab:='8_5q0']
  comps_ad1[[rd]] <- setDT(smartbind(comps_ad1[[rd]],q5))

}



## --------------------------------------------------------
## LOESS THE RAW DATA FOR A SMOOTH COMPARISON
## --------------------------------------------------------

## get all possible age and nid combos, loop through and get smoothed estimate of the raw
for(rd in rds){
  message(rd)
  combos    <- unique(comps[[rd]][,c('ab','nid'),with=FALSE])
  d         <- data.table()
  message('nat')
  for(i in 1:nrow(combos)){
    tmp <- subset(comps[[rd]], ab==combos[['ab']][i] & nid == combos[['nid']][i])
    l   <- loess(haz_dat~year, data=tmp, weight=entering, span = 0.85, family = "gaussian")
    tmp$loess_haz_dat <- predict(l,newdata=tmp)
    tmp$loess_haz_dat[tmp$loess_haz_dat<0]=0
    d <- rbind(d, tmp)
  }
  comps[[rd]] <- d
  comps[[rd]][,obs := .N,by=.(nid,ab)]
}

## --------------------------------------------------------




## --------------------------------------------------------
## CALCULATE DEVIATIONS IN ERRORS
## --------------------------------------------------------
comps0 <- copy(comps)

## ----------------------------
## FIRST FOR MORTALITY HAZARD GET ROWWISE ERROR METRICS. DEFINE AS FUNCTION
calc_rowwise_metrics <- function(data, binstoo=TRUE){

  ## error for each row
  data[, raw_err := haz_est - haz_dat]
  data[, lss_err := haz_est - loess_haz_dat]

  ## relative error for each row, raw zeros get lost on the collapse
  data[, rel_raw_err := haz_est/haz_dat]
  data[, rel_lss_err := haz_est/loess_haz_dat]
  data$rel_raw_err[data$rel_raw_err==Inf] <-NA
  data$rel_lss_err[data$rel_lss_err==Inf] <-NA

  ## absolute percentage error
  data[, pct_raw_err := abs(raw_err)/haz_dat * 100 ]
  data[, pct_lss_err := abs(lss_err)/loess_haz_dat * 100 ]
  data$pct_raw_err[data$pct_raw_err==Inf] <-NA
  data$pct_lss_err[data$pct_lss_err==Inf] <-NA

  ## ----------------------------
  ## also get errors in number entering each bin vs EEB
  if(binstoo==TRUE){
    data[, raw_EBerr     := EEB - entering] # error for each row
    data[, rel_raw_EBerr := EEB/entering]   # relative error for each row
    data[, pct_raw_EBerr := abs(raw_EBerr)/entering * 100 ] # absolute percentage error
  } else {
    data$EEB <- data$raw_EBerr <- data$rel_raw_EBerr <- data$pct_raw_EBerr <- 1
  }
  return(data)
}


## -----------------------
## COLLAPSE ERRORS TO GET SUMMARY PREDICTIVE VALIDITY METRICS
collapse_errors <- function(data, byvars){
  res <- data[,.( mean_haz   = stats::weighted.mean(loess_haz_dat,  w = entering),
                  mean_enter = mean(entering, na.rm = TRUE),
                  # mean absolute percent error
                  MdAPE_raw   = Hmisc::wtd.quantile(pct_raw_err, probs = 0.5, w = entering, na.rm = TRUE),
                  MdAPE_lss   = Hmisc::wtd.quantile(pct_lss_err, probs = 0.5, w = entering, na.rm = TRUE),
                  # mean error (bias)
                  ME_raw     = stats::weighted.mean(raw_err, w = entering),
                  ME_lss     = stats::weighted.mean(lss_err, w = entering),
                  # median relative error (bias, relative )
                  MdRE_raw    = Hmisc::wtd.quantile(rel_raw_err, probs = 0.5, weight = entering, na.rm = TRUE),
                  MdRE_lss    = Hmisc::wtd.quantile(rel_lss_err, probs = 0.5, weight = entering, na.rm = TRUE),
                  # SD of the residuals (spread of error across countries)
                  SD_res_raw = sqrt(Hmisc::wtd.var(raw_err, weight = entering)),
                  SD_res_lss = sqrt(Hmisc::wtd.var(lss_err, weight = entering)),
                  # R^2, variance explained
                  rsq_raw    = corr(cbind(haz_est,haz_dat), w = entering)^2,
                  rsq_lss    = corr(cbind(haz_est,loess_haz_dat), w = entering)^2,
                  ## now for the EB errors
                  # note here: ignore first year because some CBH children are dropped for censoring,
                  #  so its an unfair comparison if we include it
                  MdAPE_EB   = median(pct_raw_EBerr[yrsprior>0.5],   na.rm = TRUE),
                  ME_EB      = mean(raw_EBerr[yrsprior>0.5],         na.rm = TRUE),
                  MRE_EB     = median(rel_raw_EBerr[yrsprior>0.5],   na.rm = TRUE),
                  SD_res_EB  = sd(raw_EBerr[yrsprior>0.5],           na.rm = TRUE)
  ), by = c(byvars)]
}

## add metrics for each rd
comps     <- lapply(comps,     calc_rowwise_metrics)
comps_ad1 <- lapply(comps_ad1, calc_rowwise_metrics)

# how many zeros at adm1 level?
nrow(comps_ad1[[FULLMODELRD]][ab!='8_5q0'&haz_dat==0])
nrow(comps_ad1[[FULLMODELRD]][ab!='8_5q0'])

## collapse them
res_ab          <- lapply(comps,collapse_errors,byvars='ab')
res_ab_reg      <- lapply(comps,collapse_errors,byvars=c('ab','region'))
res_ab_yrsprior <- lapply(comps,collapse_errors,byvars=c('ab','yrsprior'))
res_ab_ad1          <- lapply(comps_ad1,collapse_errors,byvars='ab')
res_ab_reg_ad1      <- lapply(comps_ad1,collapse_errors,byvars=c('ab','region'))
res_ab_yrsprior_ad1 <- lapply(comps_ad1,collapse_errors,byvars=c('ab','yrsprior'))


## SAVE COLLAPSED PV INFO
for(rd in rds){
  write.csv(comps[[rd]],           sprintf('%s/%s/comp_table.csv',root,rd))
  write.csv(comps_ad1[[rd]],       sprintf('%s/%s/comp_ad1_table.csv',root,rd))

  write.csv(res_ab[[rd]],          sprintf('%s/%s/pv_table_by_ab.csv',root,rd))
  write.csv(res_ab_reg[[rd]],      sprintf('%s/%s/pv_table_by_ab_region.csv',root,rd))
  write.csv(res_ab_yrsprior[[rd]], sprintf('%s/%s/pv_table_by_ab_yrsprior.csv',root,rd))

  write.csv(res_ab_ad1[[rd]],          sprintf('%s/%s/pv_table_by_ab_ad1.csv',root,rd))
  write.csv(res_ab_reg_ad1[[rd]],      sprintf('%s/%s/pv_table_by_ab_region_ad1.csv',root,rd))
  write.csv(res_ab_yrsprior_ad1[[rd]], sprintf('%s/%s/pv_table_by_ab_yrsprior_ad1.csv',root,rd))
}
## --------------------------------------------------------





## --- FIGURES AND TABLES



## --------------------------------------------------------
## FORMAT TABLE 1
for(agg in c('natl','ad1')){
  if(agg=='natl')tab1 <- copy(res_ab[[FULLMODELRD]])
  if(agg=='ad1') tab1 <- copy(res_ab_ad1[[FULLMODELRD]])
  cols <- c('ab','mean_haz','ME_lss','SD_res_lss','MdRE_lss','MdAPE_lss','rsq_lss')
  tab1 <- tab1[, ..cols]
  setnames(tab1, cols, c('Agebin','q_a','ME','SDE','MdRE','MdAPE','Rsq'))
  write.csv(tab1,sprintf('%s/table1_%s.csv',finaloutputs,agg))
}

## --------------------------------------------------------
## MAKE FIGURE 3 -- scatter comparison

# Names of a bins clean up
abs <- data.table(ab  = c("1_nn", "2_pnn1", "3_pnn2", "4_1y", "5_2y", "6_3y", "7_4y", "8_5q0"),
                  age = as.factor(c('Neonatal (0-1 month)', 'PNN1 (1-6 months)', 'PNN2 (6-12 months)',
                                    '1-2 year (1q1)', '2-3 year (1q2)', '3-4 year (1q3)', '4-5 year (1q4)', 'All Under-5 (5q0)')))


# do it for national and ad1
for(agg in c('ad1','natl')){

  if(agg=='natl')
    comp <- comps[[FULLMODELRD]]
  if(agg=='ad1')
    comp <- comps_ad1[[FULLMODELRD]]

  comp <- merge(comp, abs, by = 'ab')
  comp$age = factor(comp$age, labels = abs$age, levels=abs$age)
  tmp <- copy(comp)

  # sort some stuff out to trick plot for nice axes
  if(agg=='natl'){

    tmp <- tmp[,.(loess_haz_dat=max(loess_haz_dat),haz_est=max(haz_est)),by=age]
    tmp[loess_haz_dat<haz_est,loess_haz_dat:=haz_est]
    tmp[loess_haz_dat>haz_est,haz_est:=loess_haz_dat]
    tmp[,haz_dat:=loess_haz_dat]
  }

  pdf(sprintf('%s/figure3_%s.pdf',finaloutputs, agg), height=6.5,width=11)
  g <-
    ggplot(comp, aes(x=loess_haz_dat, y=haz_est)) +
      geom_point(col='black',alpha=0.40,size=0.6) +
      geom_point(data=tmp,col='white',alpha=0) +  # for the scales to align
      geom_point(aes(x=0,y=0),col='white',alpha=0) +  # for the scales to align
      geom_abline(intercept=0,slope=1,color='red') +
      theme_bw() + facet_wrap(~age, scales = 'free', ncol = 4) +
      xlab('Validation Data') + ylab('Out of Sample Estimate') +
      theme_minimal()  +
      theme(strip.background = element_rect(fill="white",color='white'),
            strip.text = element_text(face="bold"),
            panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black",size=.1),
            axis.text=element_text(size=7),
            axis.ticks = element_line(size = 0.2),
            panel.border = element_blank())
  plot(g)
  dev.off()

}



## --------------------------------------------------------
## MAKE FIGURE 4 -- scatter EEB

cols<-c('#311613','#552F3E','#615673','#49849C','#37B2A2','#80D989','#EEF46D')

pdf(sprintf('%s/figure4.pdf',finaloutputs), height=6.5,width=8.5)
ggplot(comp[ab!='8_5q0'&yrsprior>0.5], aes(x=entering,y=EEB,color=age)) +
  geom_point(size=0.5) +
  geom_abline(intercept=0,slope=1,color='red') +
  theme_bw() +
  theme_minimal()  +
  scale_x_log10(limits=c(100,max(comp$entering))) +
  scale_y_log10(limits=c(100,max(comp$entering))) +
  ylab('Prediction (Expected entering age bin, EEB)') +
  xlab('Observed entering age bin in validation data') +
  scale_colour_manual(values=cols) +
  theme(strip.background = element_rect(fill="white",color='white'),
        panel.background = element_blank(), axis.line = element_line(colour = "black",size=.1),
        axis.text=element_text(size=7),
        axis.ticks = element_line(size = 0.2),
        panel.border = element_blank())
dev.off()

# this r^2 is used in the figure caption, its usually 0.96 or 0.97
summary(lm(EEB~entering,data=comp[ab!='8_5q0']))$r.squared





## --------------------------------------------------------
## MAKE FIGURE 5 and All country time series plots.

col1 <- '#072923'
col2 <- '#2C8F8F'

# nicely formatted country names instead of ISO3 codes
countries <- fread('/share/geospatial/royburst/sbh_cbh/data/locsforplots.csv')[,c('country','loc_name'),with=FALSE]
comp <- merge(comp,countries, by='country')


pdf(sprintf('%s/oos_plots_by_country.pdf',finaloutputs), height=6.5,width=11)

# loop over all testing surveys
for(n in unique(comp$nid)){
  message(n)
  # subset data
  dd <- subset(comp, nid  == n)

  # if the most recent year has too small a sample size, drop it from plotting bc its likely not representative (ie one cluster got surveyed in jan)
  tmp         <- dd[,.(N=sum(EEB)),by=.(yrsprior)]
  outlier_val <- mean(tmp[yrsprior!=min(yrsprior)]$N) - 2 *sd(tmp[yrsprior!=min(yrsprior)]$N)
  if(tmp[yrsprior==min(yrsprior)]$N < outlier_val)
    dd <- subset(dd, yrsprior > min(tmp$yrsprior))

  # get some graphing parameters
  par(mfrow=c(2,4),oma=c(0,0,2,0))
  mx1 <- max(c(dd$haz_upper[dd$age!='All Under-5 (5q0)'],dd$haz_dat[dd$age!='All Under-5 (5q0)']))*1000
  mx2 <- max(c(dd$haz_upper[dd$age=='All Under-5 (5q0)'],dd$haz_dat[dd$age=='All Under-5 (5q0)']))*1000

  # plot over each age bin
  for(a in abs$age){
    d <- subset(dd, age == a)
    plot(d$year, d$haz_est*1000, type='n', ylim=c(0, ifelse(a=='All Under-5 (5q0)',mx2,mx1)),
         ylab = 'Deaths per 1000 livebirths', xlab = 'year', main = a)
    polygon(c(d$year, rev(d$year)), c(d$haz_upper*1000, rev(d$haz_lower*1000)),
            density = 35, lwd = .5, angle = 120, col=alpha(col1,.75), border=FALSE)
    lines(d$year, d$haz_dat*1000,lwd=2, col = alpha(col2,1))
  }
  legend("bottomright",legend=c("Prediction, 95% Uncertainty Bounds","Observed Validation Data"),bty="n",
         col=c(col1,col2),density=c(35,NA),fill=c(col1,NA),border=c(NA,NA),
         lwd=c(NA,2), pt.bg = c(NA,NA), x.intersp=c(-.5,1), cex=.85)
  box(lwd=2)
  title(sprintf('Out of Sample Predictions and Validation Data for %s %s',  dd$loc_name[1],  dd$svyyr[1]), outer=TRUE)

}
dev.off()

#










## --------------------------------------------------------
## PLOTS TO COMPARE 5Q0 METHODS
## --------------------------------------------------------


## Load in the q5 data estimated from the model
q5 <- subset(comps[[FULLMODELRD]], ab == '8_5q0')

## Load in brass and IHME
brass <- fread('/share/geospatial/royburst/sbh_cbh/data/brass_q5_estimates_mab_and_tsfb_newmlt.csv')
ihm   <- fread('/share/geospatial/royburst/sbh_cbh/data/ihme_sbh_q5.csv')

# seperate out the two brass approaches
b1 <- subset(brass, method == 'mab')
b2 <- subset(brass, method == 'tsfb')

# ihme data does not have nids, so base on svyyear +/-1 and country name to match
keep <- unique(q5[,c('country','nid','svyyr'),with=FALSE])
ihm$nid <- NA
for(i in 1:nrow(keep)){
  nidyr <- keep$svyyr[i]
  cntry <- keep$country[i]
  ihm$nid[ihm$ihme_loc_id == cntry & (ihm$source_date >= (nidyr-1) & ihm$source_date <= (nidyr+1))] <- keep$nid[i]
}

# keep overlap only
nid_keep <- unique(q5$nid)[unique(q5$nid) %in% unique(ihm$nid)]

## make all data long so we can overplot them easily as time trends
q5_dat <- q5[,c('haz_dat','loess_haz_dat','year','nid','entering'),with=FALSE]
setnames(q5_dat,c('haz_dat','loess_haz_dat'),c('rawq','q'))
q5_est <- q5[,c('haz_est','haz_lower','haz_upper','year','nid'),with=FALSE]
setnames(q5_est,c('haz_est','haz_lower','haz_upper'),c('q','q_l','q_u'))
b1  <- b1[,c('q5','time','nid'),with=FALSE]
setnames(b1,c('q5','time'),c('q','year'))
b2  <- b2[,c('q5','time','nid'),with=FALSE]
setnames(b2,c('q5','time'),c('q','year'))
ihm    <- ihm[,c('q5','t','nid'),with=FALSE]
setnames(ihm,c('q5','t'),c('q','year'))
q5_dat$source <- 'Validation Data'
q5_est$source <- 'New Method'
b1$source  <- 'Standard Indirect - MAC'
b2$source  <- 'Standard Indirect - TSFB'
ihm$source    <- 'GBD-combined'
ihm[,q := q/1000]
q5 <- setDT(do.call(smartbind,list(q5_dat,q5_est,b1,b2,ihm)))
q5 <- merge(q5,keep,by='nid',all.x=TRUE)

# commented this row out below, in response to RTR, we keep all surveys
#  even if we dont have the GBD-combined estimates for them
#q5 <- subset(q5,nid %in% nid_keep)

# limit to past 15 years and 1990
dat <- subset(q5, year >= 1990)
dat[,max_year := max(year),by=nid]
dat <- subset(dat, max_year-year <= 15)
dat <- subset(dat,!is.na(country))

# svy name
dat[, svy := paste(country,svyyr)]

g1 <- ggplot(dat,aes(x=year,y=q*1000,color=source,group=source))+
  geom_line(linetype=1,lineend='round')+
  ggtitle('') + ylim(0,300) +
  ylab('5q0, Deaths per 1000') + xlab('Years') +
  facet_wrap(~svy, scales='free') +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_manual(values=c("#408F68", "#F4AC56", "#484A62", "#0061ff", "red"))

pdf(sprintf('%s/all_q5_comparisons_suppfig3.pdf',finaloutputs), height=20,width=15)
plot(g1)
dev.off()







## --
## TABLE OF PV METRICS FOR 5Q0 STUFF..

# align the dates ( could use linear interpolation instead to be a bit better here)
ihm[,   year := as.numeric(paste0(substr(as.character(year),1,4),'.','5'))]
b1[,    year := as.numeric(paste0(substr(as.character(year),1,4),'.','5'))]
b2[,    year := as.numeric(paste0(substr(as.character(year),1,4),'.','5'))]

# fixz some pesky names things
setnames(q5_dat,c('rawq','q'),c('haz_dat','loess_haz_dat'))

for(d in c('ihm','b1','b2', 'q5_est')){

  assign(d,merge(get(d),q5_dat  ,by   = c('nid','year')))
  setnames(get(d),'q','haz_est')

  assign(d,merge(get(d),keep,by='nid'))
  get(d)[, yrsprior := svyyr-year]
  #assign(d, subset(get(d), nid %in% nid_keep))

  # Use the functions from ./get_summary_pv_metrics.R
  assign(d,calc_rowwise_metrics(get(d),binstoo=FALSE))

  assign(paste0(d,'_all'),collapse_errors(get(d),byvars='source.x'))
  assign(paste0(d,'_yrsprior'),collapse_errors(get(d),byvars=c('source.x','yrsprior')))


  #
  assign(paste0(d,'_5yr'),collapse_errors(get(d)[yrsprior<5,],byvars='source.x'))
  assign(paste0(d,'_5to15yr'),collapse_errors(get(d)[yrsprior>5&yrsprior<15,],byvars='source.x'))

}

cleantab <- function(tab1){
  cols <- c('source.x','mean_haz','ME_lss','SD_res_lss','MdRE_lss','MdAPE_lss','rsq_lss')
  tab1 <- tab1[, ..cols]
  setnames(tab1, cols, c('method','q_a','ME','SDE','MdRE','MdAPE','Rsq'))
}


# different numbers of obs make it hard to interpret these, use figure instead
firstfive  <- cleantab(rbind(q5_est_5yr, ihm_5yr, b1_5yr, b2_5yr ))
lastten    <- cleantab(rbind(q5_est_5to15yr, ihm_5to15yr, b1_5to15yr, b2_5to15yr ))
allyears   <- cleantab(rbind(q5_est_all, ihm_all, b1_all, b2_all ))

write.csv(firstfive,  sprintf('%s/table2_pv_table_compare_5q0_firstfiveyears_noniddrop.csv', finaloutputs))
write.csv(lastten,    sprintf('%s/table2_pv_table_compare_5q0_lasttenyears_noniddrop.csv', finaloutputs))
write.csv(allyears,    sprintf('%s/table2_pv_table_compare_5q0_allyears_noniddrop.csv', finaloutputs))



## PLOT BY YRS PRIOR TO SURVEY
dat <- do.call(smartbind,list(ihm_yrsprior,b1_yrsprior,b2_yrsprior,q5_est_yrsprior))
setnames(dat,'source.x','model')

g1 <- ggplot(dat, aes(x=yrsprior, y=ME_lss, colour=model,group=model))+
  geom_line() + ylab('Mean Error (ME)') +
  xlab('Years Prior to Survey') + theme_bw()+
  scale_color_manual(values=c("#408F68", "#F4AC56", "#484A62", "#0061ff"))+
  theme(legend.position="none") +
  geom_abline(col='red',intercept=0,slope=0,linetype=2)
g2 <- ggplot(dat, aes(x=yrsprior, y=SD_res_lss, colour=model,group=model))+
  geom_line() + ylab('Standard Deviation of the Errors (SDE)') +
  xlab('Years Prior to Survey') + theme_bw()+
  scale_color_manual(values=c("#408F68", "#F4AC56", "#484A62", "#0061ff"))+
  theme(legend.position="none") + ylim(-0.001,0.05) +
  geom_abline(col='red',intercept=0,slope=0,linetype=2)
g3 <- ggplot(dat, aes(x=yrsprior, y=MdRE_lss, colour=model,group=model))+
  geom_line() + ylab('Median Relative Error (MRE)') +
  xlab('Years Prior to Survey') + theme_bw()+
  scale_color_manual(values=c("#408F68", "#F4AC56", "#484A62", "#0061ff"))+
  theme(legend.position="none") +
  geom_abline(col='red',intercept=1,slope=0,linetype=2)
g4 <- ggplot(dat, aes(x=yrsprior, y=MdAPE_lss, colour=model,group=model))+
  geom_line() + ylab('Median Absolute Percentage Error (MAPE)') +
  xlab('Years Prior to Survey') + theme_bw()+
  scale_color_manual(values=c("#408F68", "#F4AC56", "#484A62", "#0061ff"))+
  theme(legend.position="none") + ylim(-0.01,101) +
  geom_abline(col='red',intercept=0,slope=0,linetype=2)
g5 <- ggplot(dat, aes(x=yrsprior, y=rsq_lss, colour=model,group=model))+
  geom_line() + ylab(expression(Proportion~of~Variance~Explained~R^{2})) +
  xlab('Years Prior to Survey') + theme_bw()+
  scale_color_manual(values=c("#408F68", "#F4AC56", "#484A62", "#0061ff")) +
  geom_abline(col='red',intercept=1,slope=0,linetype=2)

pdf(sprintf('%s/q5_metrics_yrsprior_allnids.pdf',finaloutputs),height=6,width=12)
g<-grid.arrange(grobs = list(g3,g4,g5), layout_matrix = t(cbind(c(rep(1,4),rep(2,4),rep(3,7)))))
dev.off()





## --------------- COMPARISON OF METHODS FOR APPENDIX

## load pv info
comp   <- data.table()
res_ab <- data.table()
for(rd in rds){
  tmpc <- fread( sprintf('%s/%s/comp_table.csv',root,rd) )
  tmpr <- fread( sprintf('%s/%s/pv_table_by_ab.csv',root,rd) )
  tmpc$rd <- rd; tmpr$rd <- rd
  tmpc$specification <- labs[rds==rd]
  tmpr$specification <- labs[rds==rd]
  comp <- rbind(comp,tmpc); res_ab <- rbind(res_ab,tmpr)
}


# rename age bins so its pretty
res_ab[,ab := as.factor(substring(ab, 3))]

res_ab$ab <- factor(res_ab$ab, levels = c("nn",'pnn1','pnn2','1y','2y','3y','4y','5q0'))
res_ab$specification <- factor(res_ab$specification,
                               levels = labs)

# make plots
g4 <- ggplot(res_ab, aes(x=ab, y=MdAPE_lss, colour=specification, group=specification))+
  #geom_point(size=1.2,alpha=.1,position=position_dodge(width=0.35)) +
  geom_point(shape=95, size= 12) +
  theme_bw() + ylim(0,50) +
  ylab('Median Absolute Percentage Error (MAPE)') + xlab('') +
  geom_abline(col='red',intercept=0,slope=0,linetype=2) +
  theme(legend.position="none")
g1 <- ggplot(res_ab, aes(x=ab, y=ME_lss, colour=specification))+
  #geom_point(size=1.2,alpha=.1,position=position_dodge(width=0.35)) +
  geom_point(shape=95, size= 12) +
  theme_bw() + ylim(-0.005,0.005) +
  ylab('Mean Error (ME)') + xlab('') +
  geom_abline(col='red',intercept=0,slope=0,linetype=2) +
  theme(legend.position="none")
g5 <- ggplot(res_ab, aes(x=ab, y=rsq_lss, colour=specification))+
  geom_point(size= 12, shape=95) +
  theme_bw() + ylim(.5,1.01) +
  ylab(expression(Proportion~of~Variance~Explained~R^{2})) + xlab('') +
  geom_abline(col='red',intercept=1,slope=0,linetype=2)
g2 <- ggplot(res_ab, aes(x=ab, y=SD_res_lss, colour=specification))+
  geom_point(size= 12, shape=95) +
  theme_bw() + ylim(-0.001,0.03) +
  ylab('Standard Deviation of the Errors (SDE)') + xlab('') +
  geom_abline(col='red',intercept=0,slope=0,linetype=2) +
  theme(legend.position="none")
g3 <- ggplot(res_ab, aes(x=ab, y=MdRE_lss, colour=specification))+
  geom_abline(col='grey',intercept=.8,slope=0,linetype=2,size=.5) +
  geom_abline(col='grey',intercept=1.25,slope=0,linetype=2,size=.5) +
  #geom_point(size=1.2,alpha=.1) +
  geom_point(data=res_ab,  shape = 95, size= 12) +
  theme_bw() + ylim(0.8,1.3) +
  ylab('Median Relative Error (MRE)') + xlab('') +
  geom_abline(col='red',intercept=1,slope=0,linetype=2) +
  theme(legend.position="none")

pdf(sprintf('%s/spec_compare_suppfig1.pdf',finaloutputs),height=6,width=12)
g<-grid.arrange(grobs = list(g3,g4,g5), layout_matrix = t(cbind(c(rep(1,5),rep(2,5),rep(3,7)))))
dev.off()
