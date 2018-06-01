## #############################################################
## Roy Burstein
## May 2018
## This script formats the external validation citation table
##  which is in the appendix of the paper
##  Table in cludes SBH and CBH sources used but outside of the 
##  DHS sources used for the training data, which are in a 
##  seperate table. 
## ##############################################################


## set up libraries
library(data.table)


## Load data
j <- ifelse(Sys.info()['sysname']=='Linux','/home/j','J:')
# pull in the info on the nids used in the external validation, 
info  <- fread('/home/j/temp/nathenry/u5m/roy_sbh_paper/new_nids/external_validation_nids.csv')
# get citation/full title names of the surveys
names <- fread(paste0(j,'/WORK/11_geospatial/02_processed data/U5M/Global/data_list.csv'))[,c('nid','full_title'), with = FALSE]
cite  <- fread('/home/j/temp/nathenry/u5m/roy_sbh_paper/new_nids/nids_with_citations_validation.csv')[,c('nid','citation'), with = FALSE]

## Pull summary information from all the external SBH data sources
if(!file.exists(paste0(j,'/WORK/11_geospatial/02_processed data/U5M/Global/sbh_summary_info.csv'))){
  dd <- readRDS(paste0(j,'/WORK/11_geospatial/02_processed data/U5M/Global/Global_SBH_1988_2017.RDS'))
  dd <- na.omit(dd[, c('nid','ceb','ced'), with = FALSE])
  # count, ceb, ced, mothers, women
  dd <- dd[, .(ceb    = sum(ceb),
               ced    = sum(ced),
               woman  = .N,
               mother = sum(ceb>0)), 
           by = nid]
  write.csv(dd,paste0(j,'/WORK/11_geospatial/02_processed data/U5M/Global/sbh_summary_info.csv'))
  sbh_dd <- dd
} else {
  sbh_dd <- fread(paste0(j,'/WORK/11_geospatial/02_processed data/U5M/Global/sbh_summary_info.csv'))
}


## Pull CBH information as well
if(!file.exists(paste0(j,'/WORK/11_geospatial/02_processed data/U5M/Global/cbh_summary_info.csv'))){
  dd <- readRDS(paste0(j,'/WORK/11_geospatial/02_processed data/U5M/Global/Global_CBH_1988_2017.RDS'))
  dd <- na.omit(dd[, c('nid','child_alive'), with = FALSE])
  # count, ceb, ced
  dd <- dd[, .(ceb    = .N,
               ced    = sum(child_alive==0)), 
           by = nid]
  write.csv(dd,paste0(j,'/WORK/11_geospatial/02_processed data/U5M/Global/cbh_summary_info.csv'))
  cbh_dd <- dd
} else {
  cbh_dd <- fread(paste0(j,'/WORK/11_geospatial/02_processed data/U5M/Global/cbh_summary_info.csv'))
}


## add on number of children in CBH and SBH surveys
info <- merge(info, sbh_dd, by ='nid', all.x = TRUE)
info <- merge(info, cbh_dd, by ='nid', all.x = TRUE)

## make the number look pretty
info[, children := ceb.x]
info[is.na(children), children := ceb.y]
info[, children := prettyNum(children,big.mark=',')]

## clean up some  variable names
info <- info[,c('nid','loc_name','data_type','children','svyyr'),with=FALSE]
setnames(info, c('loc_name','data_type'), c('country','questionnare'))

## add survey name/citation
info <- merge(info,names,by='nid',all.x=TRUE)
info <- merge(info,cite,by='nid',all.x=TRUE)

# clean up and finalize
info[, citation := gsub('<p>', '',citation)]
info[, citation := gsub('</p>','',citation)]

# odd ones that were missed in the citaiton pull
info[nid == 237958, nid := 22125]
info[nid == 22125,  citation := "Statistics Botswana, Botswana Familiy Health Survey 2007-2008"]
info[nid == 6825, nid := 6842]
info[nid == 6842,  citation := "Statistics Indonesia, Indonesia National Socioeconomic Survey 2001"]
info[nid == 23165, nid := 23183 ]
info[nid == 23183 ,  citation := "International Institute for Population Sciences (India). India District Level Household Survey 1998-1999. Mumbai, India: International Institute for Population Sciences (India)."]
info[nid == 7140,  country := 'Jamaica']
info[nid == 39450, country := 'Jamaica']
info[nid == 7149,  country := 'Jamaica']
info[nid == 39396, country := 'Iran']

info <- info[order(country,svyyr)]
 
## finally, save the table 
write.csv(info,'/share/geospatial/royburst/sbh_cbh/plots/summary_info_data_additional_external_validation.csv')





