## #####################################################################################################
## Roy Burstein
## Replicate the standard indirect method, mother age and time since first birth variants.
## Use IGME MLT matchei
## #####################################################################################################


library(data.table)

####################################################################################################################
# 1. Clean and prepare the data and MLT inputs

##############
# load in
# load in the DHS data
full_dat <- setDT(readRDS("/home/j/WORK/11_geospatial/02_processed data/U5M/Global/CBH_with_SBH_1988_2017.RDS"))

# load in a table showing which MLT is used by IGME (sourced from www.childmortality.org)
igme_mlt <- fread('/share/geospatial/royburst/sbh_cbh/data/igme_mlt.csv')

# load in UN MLTs downloaded from https://esa.un.org/unpd/wpp/Download/Other/MLT/
m <- fread('/share/geospatial/royburst/sbh_cbh/data/MLT_UN2011_130_1y_complete.csv')


##############
# DHS clean
# Drop NAs Just three surveys with lots of NAs, not used here anyway because not recent for comparison (112500, 20351, 43016) some not even DHS
full_dat <- subset(full_dat, ! nid %in% c(112500, 20351, 43016))

# make a mother id
full_dat[, mid := paste0(nid,'-',trim(mid))]
setkey(full_dat, 'mid')

# drop mothers with missing ceb or ceb (only 3)
full_dat[, drop := sum(is.na(ceb) + is.na(ced)), by = mid] # table(full_dat$drop)
full_dat <- subset(full_dat, drop == 0)

# get a number for time since first birth for each woman
full_dat <- full_dat[, fb_cmc   := min(child_dob_cmc), by = mid]
full_dat <- full_dat[, tsfb_yr  := (interview_date_cmc - fb_cmc)/12]

# identify the tsfb grouping
full_dat[, tsfb_bin := findInterval(tsfb_yr,c(0,5,10,15,20,25))] # table(floor(full_dat$tsfb_yr),full_dat$tsfb_bin

# make it one woman per row
full_dat <- unique(full_dat)

# Scaled survey weights so they sum to the correct numbers of individuals
full_dat[, weight := weight/sum(weight,na.rm=TRUE)*(.N), by = nid]

# Get SBH summaries for this data both by mothers age bin (MAB) and time since first birth bin (TSFB)
# by mothers age
mab  <- full_dat[,.(ceb = sum(ceb*weight), ced = sum(ced*weight), num_wmn = sum(weight)),
                 by = .(nid, country, year, mothers_age_group)]
# by time since first birth
tsfb <- full_dat[,.(ceb = sum(ceb*weight), ced = sum(ced*weight), num_wmn = sum(weight)),
                 by = .(nid, country, year, tsfb_bin)]

# get ferility ratios needed for the respective methods
tmp <- copy(mab)
tmp[, parity := ceb/num_wmn]
tmp <- subset(tmp, mothers_age_group %in% c(1,2,3))[,c('nid','mothers_age_group','parity'), with=FALSE]
tmp <- reshape(tmp, idvar = "nid", timevar = "mothers_age_group", direction = "wide")
tmp[, p1p2 := parity.1/parity.2]
tmp[, p2p3 := parity.2/parity.3]

# check all is good, these should be empty
tmp[is.na(parity.1)]
tmp[is.na(parity.2)]
tmp[is.na(parity.3)]


# merge back on those fertility ratios
mab  <- merge(mab,  tmp[,c('nid','p1p2','p2p3'),with=FALSE], by='nid', all.x = TRUE)
tsfb <- merge(tsfb, tmp[,c('nid','p1p2','p2p3'),with=FALSE], by='nid', all.x = TRUE)

# combine them
mab[,  method  := 'mab']
tsfb[, method := 'tsfb']
setnames(mab, 'mothers_age_group', 'bin')
mab[, bin := as.numeric(bin)]
setnames(tsfb, 'tsfb_bin', 'bin')
mab <- subset(mab, bin %in% 1:7)
tsfb <- subset(tsfb, bin %in% 1:5)
sd <- rbind(mab, tsfb)

# cd.ceb ratio, start using Manual X notation
sd[, Di := ced/ceb]



##############
# clean the information on MLTs used by IGME to try and match as best we can to what they do
# for those non standard ones or B3, use the west mlt
igme_mlt[! family %in% c('East','West','North','South','UN General','UN Latin American'), family := 'West']

# seperate column for coeffs (keep to CD to match whats available in Manual X (this only affects two countries with the UN MLTs))
igme_mlt[, family_coef := tolower(family)]
igme_mlt[! family_coef %in% c('east','west','north','south'), family_coef := 'west']


# merge on the family info
sd <- merge(sd, igme_mlt, by = 'country')



##############
# Clean up the UN MLTs in a way that is usable
# keep only needed variables and ages
m <- m[, c('Family','Sex','E0','age','qx1'), with = FALSE]
m <- subset(m, age <= 20)

# combine sexes (simple mean for both)
m <- m[, .(q = mean(qx1)), by = .(Family,E0,age)]

# make q cumulative (ie Xq0)
a <- m$age
m <- m[, .(X = 1 - cumprod(1 - q)), by = .(Family,E0)]
m[, age := a]

# reshape it like a matrix for easy search
qx_table <- reshape(m, timevar = 'E0', idvar = c('age','Family'), direction = 'wide')
setnames(qx_table, 'Family', 'family')

# add upper limit
qx_table$X101 <- 0

# fix some names so they match data
qx_table[family == 'Latin',   family := 'UN Latin American']
qx_table[family == 'General', family := 'UN General']

# save vector of data columns to use later in the search
dat_cols      <-  grep('X',colnames(qx_table))

# save prepped mlt stuff for later
write.csv(igme_mlt,'/share/geospatial/royburst/sbh_cbh/data/prepped_igme_mlt.csv')
write.csv(qx_table,'/share/geospatial/royburst/sbh_cbh/data/prepped_MLT_UN2011_130_1y_complete.csv')



####################################################################################################################
# 2. Calculate indirect

# Get coefficient tables

# mother age coefficents from Manual X
a_cm_coeffs <- setDT(readRDS('/share/geospatial/royburst/sbh_cbh/data/manualX_table47.RDS'))
a_rt_coeffs <- setDT(readRDS('/share/geospatial/royburst/sbh_cbh/data/manualX_table48.RDS'))

# tsfb coefficients from Moultrie Chapter 16
b_cm_coeffs <- setDT(fread('/share/geospatial/royburst/sbh_cbh/data/moultrie_table16.3.csv'))
b_rt_coeffs <- setDT(fread('/share/geospatial/royburst/sbh_cbh/data/moultrie_table16.4.csv'))


# Equation: k(i) = a(i) + b(i)(P(1)/P(2)) + c(i)(P(2)/P(3))
tmpa <- merge(subset(sd,method=='mab'),
              a_cm_coeffs[,c('family','index','coef_a','coef_b','coef_c'),with=FALSE],
              by.x = c('family_coef','bin'), by.y = c('family','index'), all.x = TRUE)
tmpb <- merge(subset(sd,method=='tsfb'),
              b_cm_coeffs[,c('family','index','coef_a','coef_b','coef_c'),with=FALSE],
              by.x = c('family_coef','bin'), by.y = c('family','index'), all.x = TRUE)

tmpa[, ki :=  coef_a + coef_b*p1p2 + coef_c*p2p3]
tmpb[, ki :=  coef_a + coef_b*p1p2 + coef_c*p2p3]


# calculate q(x) for each child age, k(i) * D(i)
tmpa[, qx := ki * Di]
tmpb[, qx := ki * Di]


# mark which age each qx is associated with
ca <- data.table(mothers_age_group = 1:7, x = c(1,2,3,5,10,15,20))
tmpa$x <- ca$x[match(tmpa$bin,ca$mothers_age_group)]

ca <- data.table(tsfb_group = 1:5, x = c(2,5,5,5,10))
tmpb$x <- ca$x[match(tmpb$bin,ca$tsfb_group)]


# for each row, we need to find the column most closely associated with it, by age (x), family, and qx
for(tmpd in c('tmpa','tmpb')){
  d <- get(tmpd)
  d$q5 <- as.numeric(NA)
  for (i in 1:nrow(d)) {
    tmp_qx      <- d[["qx"]][i]
    if(tmp_qx <= 0){
      q5_est = 0
    } else {
      # identify the columns
      tmp_mlt     <- subset(qx_table, family == d[['family']][i] & age == d[['x']][i])
      tmp_mlt_q5  <- subset(qx_table, family == d[['family']][i] & age == 5)
      poss_vals   <- as.numeric(as.matrix(tmp_mlt)[,dat_cols])
      names(poss_vals) <- dat_cols
      closest_val <- which(abs(poss_vals-tmp_qx)==min(abs(tmp_qx-poss_vals)))
      closest_higher <- poss_vals[closest_val] > tmp_qx
      col_hi <- as.numeric(ifelse(closest_higher,names(closest_val),as.numeric(names(closest_val))-1))
      col_lo <- col_hi+1

      # sort out the q5 tranformation
      qx_lo <- as.numeric(as.matrix(tmp_mlt)[,col_lo])
      qx_hi <- as.numeric(as.matrix(tmp_mlt)[,col_hi])
      q5_lo <- as.numeric(as.matrix(tmp_mlt_q5)[,col_lo])
      q5_hi <- as.numeric(as.matrix(tmp_mlt_q5)[,col_hi])
      q5_est <- q5_lo + (tmp_qx - qx_lo) / (qx_hi - qx_lo) * (q5_hi - q5_lo)
    }
    d$q5[i] <- q5_est
  }
  assign(tmpd,d)
}


# calculate t(x) time adjustment for each maternal age group, based on table 48
# Equation: t(x) = a(i) + b(i)(P(1)/P(2)) + c(i)(P(2)/P(3))
tmpa <- merge(tmpa,
              a_rt_coeffs[,c('family','index','coef_a2','coef_b2','coef_c2'),with=FALSE],
              by.x = c('family_coef','bin'), by.y = c('family','index'), all.x = TRUE)
tmpb <- merge(tmpb,
              b_rt_coeffs[,c('family','index','coef_a2','coef_b2','coef_c2'),with=FALSE],
              by.x = c('family_coef','bin'), by.y = c('family','index'), all.x = TRUE)

tmpa[, tx := coef_a2 + coef_b2*p1p2 + coef_c2*p2p3]
tmpa[, time := year - tx]

tmpb[, tx := coef_a2 + coef_b2*p1p2 + coef_c2*p2p3]
tmpb[, time := year - tx]


# combine them and save
final <- rbind(tmpa,tmpb)

# save the output
final <- final[order(nid,country,time,method)]
write.csv(final,file = '/share/geospatial/royburst/sbh_cbh/data/brass_q5_estimates_mab_and_tsfb_newmlt.csv')



# do a plot comparison of the two methods
pdf('/share/geospatial/royburst/sbh_cbh/data/brass_methods_compare_newmlt.pdf',height=8,width=10)
for(n in unique(final$nid)){
  message(n)
  tmp <- subset(final,nid==n)
  title <- paste(tmp$country[1], tmp$year[1], tmp$nid[1])
  g <- ggplot(tmp, aes(x=time,y=qx,group=method,color=method)) +
    theme_bw() + geom_line() + ggtitle(title)
  plot(g)
}
dev.off()
