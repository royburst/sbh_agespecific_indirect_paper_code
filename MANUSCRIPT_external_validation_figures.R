## Roy Burstein
## 2018
## Plots and tables to do with external validation


## ---------------------------------------------------------------------------
## Do some basic setup stuff

# Load packages
library(data.table)
library(ggplot2)
library(gtools)
library(ggrepel)
library(gridExtra)
library(Hmisc)
library(dplyr)
library(parallel)

# change root of data share drive depending on which OS using.
j <- ifelse(Sys.info()['sysname']=='Linux','/home/j','J:')

# where to save plots
FULLMODELRD  <-  '2018_03_08_09_56_20'
finaloutputs <- paste0('/share/geospatial/royburst/sbh_cbh/plots/',FULLMODELRD)



## ---------------------------------------------------------------------------
## load data
sbh <- readRDS(paste0(j,'/temp/royburst/sbh_cbh/verification_data/prepared_verification_sbh_tmp.RDS'))
cbh <- readRDS(paste0(j,'/temp/royburst/sbh_cbh/verification_data/prepared_verification_cbh_1988_2017_censfix.RDS'))



## ---------------------------------------------------------------------------
## Some slight cleaning is needed for the cbh dataset
## collapse N and died to the survey year

# The india 2000 DHS has been split into subnational NIDS, just rename the source id here
cbh <- cbh[country=='IND'&svyyr==2000&source=='MACRO_DHS',nid:=19950]


cbh <- cbh[,.(died = sum(died), N = sum(N)),
           by = .(country,nid,source,svyyr,period,year,ab,data_type)]

# estimate hazard
cbh[, haz := died/N]

# Combine age bins to get estimates for NN, INF, 5q0
inf <- subset(cbh, ab %in% c('NN','PNN1','PNN2'))[, .(haz = 1 - prod(1-haz), N = max(N)),
                                                  by =  .(country,nid,source,svyyr,period,year,data_type)]
u5  <- cbh[,.(haz = 1 - prod(1-haz), N = max(N)),
           by =  .(country,nid,source,svyyr,period,year,data_type)]
nn  <- subset(cbh, ab %in% 'NN')
nn[,  ab := 'NN']; nn[, died := NULL]
inf[, ab := '1q0']
u5[,  ab := '5q0']
cbh <- do.call('rbind', list(nn,inf,u5))

# we'll start in 1990
cbh <- subset(cbh, year >= 1990)

# combine sbh and cbh
d <- setDT(smartbind(sbh,cbh))
d[, source := paste0(source,'_',svyyr,'_',nid)]

d[ab == 'NN',  ab := '1_NN']
d[ab == '1q0', ab := '2_1q0']
d[ab == '5q0', ab := '3_5q0']

# helps for plotting
d$f_data_type <- factor(d$data_type)

## ---------------------------------------------------------------------------
# how many ceb, ced, mothers from each of these surveys?
summ  <-  fread(paste0(j,'/WORK/11_geospatial/02_processed data/U5M/Global/sbh_summary_info.csv'))
summ  <-  subset(summ, nid %in% unique(d[data_type=='sbh']$nid))
sum(summ$ceb)     # ~153M
sum(summ$ced)     # ~11M
sum(summ$mother)  # ~54M



## ---------------------------------------------------------------------------
# loop over countries and plot data from each one

# Pull in additional info to get nice survey names for plotting
info <- fread(paste0(j,'/WORK/11_geospatial/02_processed data/U5M/Global/data_list.csv'))[,c('nid','full_title'), with = FALSE]
d <- merge(d,info,by='nid',all.x=TRUE)

# clean up a few things
d$full_title <- gsub("Multiple Indicator Cluster Survey", "MICS", d$full_title )
d$full_title <- gsub("Demographic and Health Survey", "DHS", d$full_title )
d$full_title <- gsub("Household", "HH", d$full_title )
d$full_title <- gsub("Population", "Pop.", d$full_title )
d$full_title <- gsub("Census", "Cens.", d$full_title )

# Limit to 15 year trends
d <- subset(d, svyyr - year < 16)


# clean things up to plot
d[, age := factor(substr(ab,3,1000), levels = c('NN','1q0','5q0'), labels = c('NN','1q0','5q0'))]
countries <- fread('/share/geospatial/royburst/sbh_cbh/data/locsforplots.csv')[,c('country','loc_name'),with=FALSE]
d <- merge(d,countries, by='country', all.x=TRUE)

# save info on all of these surveys to be used in comparison testing
extsvyinfo <- unique(d[,c('country','loc_name','nid','svyyr','data_type'),with=FALSE])
inf <- fread('/share/geospatial/royburst/sbh_cbh/plots/summary_info_data.csv') # DHS surveys
extsvyinfo <- extsvyinfo[! nid %in% inf$ghdx_nid,]
write.csv(extsvyinfo,'/home/j/temp/nathenry/u5m/roy_sbh_paper/new_nids/external_validation_nids.csv')

# do them plots
cnts <- sort(unique(d$loc_name))
pdf(paste0(finaloutputs,'/suppfig_SBH_external_verification_by_country.pdf'),width=12,height=8)
for(c in cnts){
  message(c)

  # subset out data by country
  tmp <- subset(d, loc_name == c)

  # make sure labels will do into tthe correct year
  tmp2 <- tmp[,mxy := max(year),by=nid]
  tmp2 <- tmp2[year == mxy]

  g <-
  ggplot(tmp, aes(x=year,y=haz*1000,group=source,color=f_data_type,fill=f_data_type)) +
    theme_bw() +
    facet_wrap(~age) +
    geom_line(size=1.5,alpha=0.5) +
    scale_colour_manual(drop   = TRUE,
                          name   = 'Data Type',
                          values = c('#F16B6F','#30A9DE'),
                          labels = c('CBH - Direct','SBH - Indirect, new method'),
                          limits = levels(tmp$f_data_type)) +
    theme(strip.background = element_blank(),
          panel.border     = element_rect(colour = "black")) +
    geom_text_repel(data=tmp2,
                    aes(label=full_title),
                    force         = 5,
                    size          = 2,
                    direction     = 'y',
                    nudge_x       = 2020,
                    box.padding   = 0.5,
                    segment.size  = 0.1,
                    segment.color = "grey30")+
    scale_x_continuous(breaks=c(1990,1995,2000,2005,2010,2015),limits=c(1988,2030)) +
    ylab('Estimated deaths per 1000 livebirths') + ylim(0,300) +
    ggtitle(tmp$loc_name[1])
  plot(g)
}
dev.off()


# howmany countries
length(unique(d[data_type=='sbh']$country))



## ---------------------------------------------------------------------------
# Do scatter plots for each SBH-CBH country-year pair.

# Format data for easy comparison
dd <- d[,c('ab','nid','year','country','haz','data_type','N'),with=FALSE]
dd_cbh <- subset(dd, data_type=='cbh')
setnames(dd_cbh, c('haz','nid'),c('haz_cbh','nid_cbh'))
dd_cbh[,data_type:=NULL]
dd_sbh <- subset(dd, data_type=='sbh')
dd_sbh[, N := NULL]
setnames(dd_sbh, c('haz','nid'),c('haz_sbh','nid_sbh'))
dd_sbh[,data_type:=NULL]

# loess cbh, as we did in main validation
combos    <- unique(dd_cbh[,c('ab','nid_cbh'),with=FALSE])
ddd       <- data.table()
for(i in 1:nrow(combos)){
  tmp <- subset(dd_cbh, ab==combos[['ab']][i] & nid_cbh == combos[['nid_cbh']][i])
  if(nrow(tmp)>2){
    l   <- loess(haz_cbh~year, data=tmp, span = 0.85, family = "gaussian")
    tmp$haz_cbh_loess <- predict(l,newdata=tmp)
    tmp$haz_cbh_loess[tmp$haz_cbh_loess<0] <- 0
  } else {
    tmp$haz_cbh_loess <- tmp$haz_cbh
  }
  ddd <- rbind(ddd, tmp)
}
dd_cbh <- ddd


# combine sbh cbh
dd <- dd_sbh %>%  left_join(dd_cbh, by = c('ab','year','country'))
dd <- setDT(dd)

# combine cbh cbh
dd2 <- dd_cbh %>%  left_join(dd_cbh, by = c('ab','year','country'))
dd2 <- setDT(dd2)
dd2 <- subset(dd2, nid_cbh.x != nid_cbh.y)
dd2[,tmp:= paste0(ab,year*round(haz_cbh.x,3)*round(haz_cbh.y,3)*nid_cbh.x*nid_cbh.y)] # also drop duplicates that were in reverse order
dd2 <- dd2[!duplicated(dd2$tmp),]

# some other figures to pull:
# country-years of data with SBH only
sum(!is.na(dd$haz_sbh)    & is.na(dd$haz_cbh)) # sbh in cy without cbh
nrow(unique(dd[!is.na(haz_sbh) & is.na(haz_cbh), c('country','year'), with=FALSE])) # number of cy without cbh but with sbh
sum(!is.na(dd$haz_sbh)    & !is.na(dd$haz_cbh))
sum(!is.na(dd2$haz_cbh.y) & !is.na(dd2$haz_cbh.y))

# clean up for plotting
dd[,  age := factor(substr(ab,3,1000), levels = c('NN','1q0','5q0'), labels = c('NN','1q0','5q0'))]
dd2[, age := factor(substr(ab,3,1000), levels = c('NN','1q0','5q0'), labels = c('NN','1q0','5q0'))]

# trick to make sure plots are square
tmp <- dd[,.(haz_cbh_loess=max(haz_cbh_loess,na.rm=T),haz_sbh=max(haz_sbh,na.rm=T)),by=age]
tmp[haz_cbh_loess<haz_sbh,haz_cbh_loess:=haz_sbh]
tmp[haz_cbh_loess>haz_sbh,haz_sbh:=haz_cbh_loess]
tmp[,haz_sbh:=haz_cbh_loess]
tmp[,haz_cbh_loess.y:=haz_cbh_loess]
tmp[,haz_cbh_loess.x:=haz_cbh_loess]

# do the plotting
g4 <- ggplot(dd, aes(x=haz_cbh_loess,y=haz_sbh)) + theme_bw() +
       geom_point(size=.4) + geom_abline(intercept=0,slope=1,color='red') +
       facet_wrap(~age,scales='free') + ggtitle('CBH versus SBH concurrent') +
       ylab('SBH') + xlab('CBH') +
       geom_point(data=tmp,col='white',alpha=0) +  # for the scales to align
       theme(strip.background = element_blank(),
            panel.border     = element_rect(colour = "black"))
g5 <-  ggplot(dd2, aes(x=haz_cbh_loess.x,y=haz_cbh_loess.y)) + theme_bw() +
       geom_point(size=.4) + geom_abline(intercept=0,slope=1,color='red') +
       facet_wrap(~age,scales='free') + ggtitle('CBH versus CBH concurrent')+
       ylab('CBH 1') + xlab('CBH 2') +
       geom_point(data=tmp,col='white',alpha=0) +  # for the scales to align
       theme(strip.background = element_blank(),
             panel.border     = element_rect(colour = "black"))

pdf(paste0(finaloutputs,'/figure7_CBH_SBH_external_scatter.pdf'),width=12,height=8)
grid.arrange(g4,g5)
dev.off()


# get R-squared in both, also while correcting for nid?
## NN, without NID correction
summary(lm(haz_sbh~haz_cbh_loess,data=subset(dd,ab=='1_NN')))$r.squared
summary(lm(haz_cbh_loess.y~haz_cbh_loess.x,data=subset(dd2,ab=='1_NN')))$r.squared

## NN, with NID correction
summary(lm(haz_sbh~haz_cbh_loess+factor(nid_sbh),data=subset(dd,ab=='1_NN')))$r.squared
summary(lm(haz_cbh_loess.y~haz_cbh_loess.x+factor(nid_cbh.y),data=subset(dd2,ab=='1_NN')))$r.squared

## 1q0, without NID correction
summary(lm(haz_sbh~haz_cbh_loess,data=subset(dd,ab=='2_1q0')))$r.squared
summary(lm(haz_cbh_loess.y~haz_cbh_loess.x,data=subset(dd2,ab=='2_1q0')))$r.squared

## 1q0, with NID correction
summary(lm(haz_sbh~haz_cbh_loess+factor(nid_sbh),data=subset(dd,ab=='2_1q0')))$r.squared
summary(lm(haz_cbh.y~haz_cbh_loess.x+factor(nid_cbh.y),data=subset(dd2,ab=='2_1q0')))$r.squared

# 5q0, without NID correction
summary(lm(haz_sbh~haz_cbh_loess,data=subset(dd,ab=='3_5q0')))$r.squared
summary(lm(haz_cbh_loess.y~haz_cbh_loess.x,data=subset(dd2,ab=='3_5q0')))$r.squared

# 5q0, with NID correction
summary(lm(haz_sbh~haz_cbh_loess+factor(nid_sbh),data=subset(dd,ab=='3_5q0')))$r.squared
summary(lm(haz_cbh_loess.y~haz_cbh_loess.x+factor(nid_cbh.y),data=subset(dd2,ab=='3_5q0')))$r.squared

## ---------------------------------------------------------------------------
# basic set of pvmetrics  comparisons from paper

dd <- na.omit(dd)

calc_rowwise_metrics <- function(data){
  ## error for each row
  data[, raw_err := haz_sbh - haz_cbh]
  data[, lss_err := haz_sbh - haz_cbh_loess]
  ## relative error for each row, raw zeros get lost on the collapse
  data[, rel_raw_err := haz_sbh/haz_cbh]
  data[, rel_lss_err := haz_sbh/haz_cbh_loess]
  data$rel_raw_err[data$rel_raw_err==Inf] <-NA
  data$rel_lss_err[data$rel_lss_err==Inf] <-NA
  ## absolute percentage error
  data[, pct_raw_err := abs(raw_err)/haz_cbh * 100 ]
  data[, pct_lss_err := abs(lss_err)/haz_cbh_loess * 100 ]
  data$pct_raw_err[data$pct_raw_err==Inf] <-NA
  data$pct_lss_err[data$pct_lss_err==Inf] <-NA
  return(data)
}

collapse_errors <- function(data, byvars){
  res <- data[,.( mean_haz_cbh   = stats::weighted.mean(haz_cbh_loess,  w = N),
                  mean_haz_sbh   = stats::weighted.mean(haz_sbh,  w = N),
                  # mean absolute percent error
                  MdAPE_raw   = Hmisc::wtd.quantile(pct_raw_err, probs = 0.5, w = N, na.rm = TRUE),
                  MdAPE_lss   = Hmisc::wtd.quantile(pct_lss_err, probs = 0.5, w = N, na.rm = TRUE),
                  # mean error (bias)
                  ME_raw     = stats::weighted.mean(raw_err, w = N),
                  ME_lss     = stats::weighted.mean(lss_err, w = N),
                  # median relative error (bias, relative )
                  MdRE_raw    = Hmisc::wtd.quantile(rel_raw_err, probs = 0.5, weight = N, na.rm = TRUE),
                  MdRE_lss    = Hmisc::wtd.quantile(rel_lss_err, probs = 0.5, weight = N, na.rm = TRUE),
                  # SD of the residuals (spread of error across countries)
                  SD_res_raw = sqrt(Hmisc::wtd.var(raw_err, weight = N)),
                  SD_res_lss = sqrt(Hmisc::wtd.var(lss_err, weight = N))
  ), by = c(byvars)]
}

x <- collapse_errors(calc_rowwise_metrics(dd), byvars='ab')
write.csv(x,paste0(finaloutputs,'/table2_external_sbh_pv_metrics.csv'))

# trick it to get CBH CBH also
setnames(dd2, c('haz_cbh_loess.x','haz_cbh_loess.y','haz_cbh.y','N.x'), c('haz_sbh', 'haz_cbh_loess','haz_cbh','N'))
x2 <- collapse_errors(calc_rowwise_metrics(dd2), byvars='ab')







## ---------------------------------------------------------------------------
## Compare ratios NN/5q0, 1q0/5q0

dd <- dd[year>1992] # because some age binsstarted staggered
rsbh <- cbind(dd[ab=='1_NN',c('year','nid_sbh','country','N','haz_sbh'),with=FALSE],
            dd[ab=='2_1q0',c('haz_sbh'),with=FALSE],
            dd[ab=='3_5q0',c('haz_sbh'),with=FALSE])
names(rsbh) <- c('year','nid_sbh','country','N_sbh','nn_sbh','q1_sbh','q5_sbh')
rsbh[, ratio_sbh_nnq5 := nn_sbh/q5_sbh]
rsbh[, ratio_sbh_q1q5 := q1_sbh/q5_sbh]
rsbh[, ratio_sbh_nnq1 := nn_sbh/q1_sbh]

rsbh <- unique(rsbh)

rcbh <- cbind(dd[ab=='1_NN',c('year','nid_cbh','country','N','haz_cbh_loess'),with=FALSE],
              dd[ab=='2_1q0',c('haz_cbh_loess'),with=FALSE],
              dd[ab=='3_5q0',c('haz_cbh_loess'),with=FALSE])
names(rcbh) <- c('year','nid_cbh','country','N_cbh','nn_cbh','q1_cbh','q5_cbh')
rcbh[, ratio_cbh_nnq5 := nn_cbh/q5_cbh]
rcbh[, ratio_cbh_q1q5 := q1_cbh/q5_cbh]
rcbh[, ratio_cbh_nnq1 := nn_cbh/q1_cbh]

rcbh <- unique(rcbh)

rat <- rcbh %>%  left_join(rsbh, by = c('year','country'))
rat <- setDT(rat)

# plot with ggplot for the 3 ratios
rs1 <- round(summary(lm(ratio_cbh_nnq5~ratio_sbh_nnq5,data=rat))$r.squared,3)
rs2 <- round(summary(lm(ratio_cbh_q1q5~ratio_sbh_q1q5,data=rat))$r.squared,3)
rs3 <- round(summary(lm(ratio_cbh_nnq1~ratio_sbh_nnq1,data=rat))$r.squared,3)

g1 <-
ggplot(rat, aes(x=ratio_cbh_nnq5,y=ratio_sbh_nnq5)) +
  theme_bw() + geom_point(size=1,alpha=0.5) +
  geom_abline(slope=1,intercept=0,color='red',size=1.15) +
  ylab('SBH') +
  xlab('CBH') +
  xlim(0,1) + ggtitle('Ratio of NN:5q0') +
  annotate("text", x = 0.08, y = max(rat$ratio_sbh_nnq5), label = paste0('R-squared = ',rs1))
g2 <-
ggplot(rat, aes(x=ratio_cbh_q1q5,y=ratio_sbh_q1q5)) +
  theme_bw() + geom_point(size=1,alpha=0.5) +
  geom_abline(slope=1,intercept=0,color='red',size=1.15)  +
  ylab('SBH') +
  xlab('CBH') +
  xlim(0,1) + ggtitle('Ratio of 1q0:5q0') +
  annotate("text", x = 0.08, y = max(rat$ratio_sbh_q1q5), label = paste0('R-squared = ',rs2))
g3 <-
ggplot(rat, aes(x=ratio_cbh_nnq1,y=ratio_sbh_nnq1)) +
  theme_bw() + geom_point(size=1,alpha=0.5) +
  geom_abline(slope=1,intercept=0,color='red',size=1.15) +
  ylab('SBH') +
  xlab('CBH')  +
  xlim(0,1) + ggtitle('Ratio of NN:1q0') +
  annotate("text", x = 0.08, y = max(rat$ratio_sbh_nnq1), label = paste0('R-squared = ',rs3))

pdf(paste0(finaloutputs,'/suppfig_cbh_sbh_ratio_comparisons.pdf'),height=12,width=8)
grid.arrange(g1,g2,g3,ncol=1,nrow=3)
dev.off()
