###Fit comparator models
##Setup
#Combine datasets
alldata2 <- ACexpmat2 <- data.frame()
for (i in ctrys_fit) {
  load(paste0("data/",i,"_rates.RData"))
  alldata2 <- rbind(alldata2,data.frame(ctry=i,alldata))
  ACexpmat2 <- rbind(ACexpmat2,data.frame(ctry=i,ACexpmat))
}
alldata2 %<>% mutate(Year=Cohort+Age) #add year column

#Construct minimum cohort data frame
ctry_joy_df <- data.frame(
  ctry=substr(joy_sel_lst %>% unlist() %>% names(),1,3),
  joy=joy_sel_lst %>% unlist(),
  red=red_sel_lst %>% unlist())
ctry_joy_df %<>% mutate(mincoh=ifelse(red==T,joy-max(agerange)-5,joy-max(agerange)-10))

#Compute CFR
#Up to jump-off year
obscfr <- data.frame(joy=rep(joys,each=nrow(alldata2)),
                     alldata2[rep(1:nrow(alldata2),length(joys)),]) %>% filter(par==0) %>%
  left_join(ctry_joy_df) %>%
  filter(Cohort>=mincoh) %>%
  group_by(ctry,Cohort,joy) %>%
  mutate(mina=min(Age),maxa=max(Age)) %>%
  filter(mina==min(agerange),maxa==max(agerange)) %>%
  filter(Year<=joy) %>%
  group_by(ctry,Cohort,joy) %>% summarise(jx=sum(x))

#Complete, including the years after the jump-off year
obscfr2 <- alldata2 %>% filter(par==0, ctry %in% ctrys_fit) %>%
  group_by(ctry,Cohort) %>%
  mutate(mina=min(Age),maxa=max(Age)) %>%
  filter(mina==min(agerange),maxa==max(agerange)) %>%
  group_by(ctry,Cohort) %>% summarise(x=sum(x))

#Create JOY-country-cohort-age empty data frame
jccadat <- data.frame(expand.grid(Age=agerange,Cohort=0:39,ctry=ctrys_fit,joy=joys)) %>%
  mutate(Cohort=Cohort+joy-max(agerange)-10,Year=Age+Cohort) %>%
  left_join(ctry_joy_df) %>% filter(Cohort>=mincoh)

#Create JOY-country-cohort empty data frame
jccdat <- data.frame(expand.grid(Cohort=0:39,ctry=ctrys_fit,joy=joys)) %>%
  mutate(Cohort=Cohort+joy-max(agerange)-10) %>%
  left_join(ctry_joy_df) %>% filter(Cohort>=mincoh)


##Freeze rate
#ASFRs
fradat <- jccadat %>% left_join(alldata2 %>% filter(par==0) %>% select(-par,-N,-n)) %>%
  group_by(joy,ctry,Age) %>%
  mutate(pred=ifelse(Year<joy,x,x[which(Year==joy)])) %>%
  mutate(sd=NA,lwr95=NA,upr95=NA,lwr90=NA,upr90=NA,lwr50=NA,upr50=NA,crps=NA,logs=NA)

#Calculate CRPS (reduces to MAE)
ind <- which(fradat$Year>fradat$joy & !is.na(fradat$x))
fradat$crps[ind] <- abs(fradat$x[ind]-fradat$pred[ind])

#CFRs
frcdat <- fradat %>% group_by(ctry,Cohort,joy) %>%
  summarise(pred=sum(pred)) %>% left_join(obscfr2) %>%
  mutate(sd=NA,lwr95=NA,upr95=NA,lwr90=NA,upr90=NA,lwr50=NA,upr50=NA,crps=NA,logs=NA)

#Calculate CRPS (reduces to MAE)
ind <- which(frcdat$Cohort+max(agerange)>frcdat$joy & !is.na(frcdat$x))
frcdat$crps[ind] <- abs(frcdat$x[ind]-frcdat$pred[ind])


##Myrskylä et al.
#ASFRs
myrskadat <- jccadat %>% left_join(alldata2 %>% filter(par==0) %>% select(-par,-N,-n)) %>%
  group_by(joy,ctry,Age) %>%
  mutate(roc=(x[which(Year==joy)]-x[which(Year==(joy-4))])/4) %>%
  mutate(pred=ifelse(Year<joy,x,x[which(Year==joy)])) %>%
  mutate(addon=roc*c(rep(0,length(Year[1]:joy[1])),1:5,rep(5,45))[1:n()]) %>%
  mutate(pred=pred+addon) %>%
  mutate(pred=ifelse(pred<0|is.na(pred),0,pred)) %>%
  select(-c("roc","addon")) %>%
  mutate(sd=NA,lwr95=NA,upr95=NA,lwr90=NA,upr90=NA,lwr50=NA,upr50=NA,crps=NA,logs=NA)

#Calculate CRPS (reduces to MAE)
ind <- which(myrskadat$Year>myrskadat$joy & !is.na(myrskadat$x))
myrskadat$crps[ind] <- abs(myrskadat$x[ind]-myrskadat$pred[ind])

#CFRs
myrskcdat <- myrskadat %>% group_by(ctry,Cohort,joy) %>%
  summarise(pred=sum(pred)) %>% left_join(obscfr2) %>%
  mutate(sd=NA,lwr95=NA,upr95=NA,lwr90=NA,upr90=NA,lwr50=NA,upr50=NA,crps=NA,logs=NA)

#Calculate CRPS (reduces to MAE)
ind <- which(myrskcdat$Cohort+max(agerange)>myrskcdat$joy & !is.na(myrskcdat$x))
myrskcdat$crps[ind] <- abs(myrskcdat$x[ind]-myrskcdat$pred[ind])


##Schmertmann et al.
#The code to fit the Schmertmann et al. (2014) model is not provided but can be obtained straightforwardly
#by the user by downloading the R code provided on the project website ("https://schmert.net/cohort-fertility/").
#Note the following minor modifications to the code:
#1. Instead of using the "VV" HFD files (period fertility rates by calendar year, age reached during
#the year and birth cohort), to be consistent with our age-cohort setup we use the "VH" files (fertility
#rates by birth cohort and age in completed years). It is necessary to update the list of country full
#names ("full.name") in "input rate and exposure data.R" to include the HFD countries that have been
#added since the Schmertmann et al. (2014) code was published. We also remove countries from the additional
#data that now have HFD data available.
#2. To fit the model for each jump-off year, we set 'forecast.year' to be the jump-off year of interest,
#thereby fitting to the appropriate (full) cohort range automatically. Consistent with this, we adjust the
#historical cohort range depending on the earliest cohort included such that if the earliest cohort is
#cohort c, the historical cohort range is set as 1900 to (c-7). If the particular country-JOY combination
#has a reduced cohort setup (red=T), we adjust the range of cohorts fitted to, as well as the historical
#cohort range, appropriately. We create two versions of each of the files, editing the reduced versions to 
#save any outputs in a file with suffix "_red", and edit all lines relating to the full cohort setup (40
#cohorts) so that they refer to the reduced 35 cohorts only.
#3. For each country-JOY combination, we run the model twice, once as normal and once with a slight change to
#the matrix "C" such that its entries corresponding to the observed age-cohort combinations (where "vis" is
#equal to 1) are set to zero and therefore not included in the CFR calculation. In this way, the CFR posterior
#summaries refer to what is added on to the rates that were observed in the fitting period, assuming that these
#rates were observed without error (we also do this for the proposed model - see lines 324-355 of
#"r/0_functions.r"). For each jump-off year and reduced/full cohort range, we save the two model outputs in the
#'output' subfolder as 'schmert_<joy>_1.dput' and 'schmert_<joy>_2.dput' for the full cohort setup, and 
#'schmert_<joy>_1_red.dput' and 'schmert_<joy>_2_red.dput' for the reduced cohort setup. There are thus a maximum
#of 4 output files for each JOY, depending on whether there are any countries with the reduced cohort setup.

#Process output
schmertout1 <- schmertout2 <- schmertout1_red <- schmertout2_red <- list()
for (i in 1:length(joys)) {
  schmertout1[[i]]     <- dget(paste0("output/schmert_",joys[i],"_1.dput"))
  schmertout2[[i]]     <- dget(paste0("output/schmert_",joys[i],"_2.dput"))
  if (ctry_joy_df %>% filter(joy==joys[i]) %>% pull(mincoh) %>% unique() %>% length() > 1)   {
    schmertout1_red[[i]] <- dget(paste0("output/schmert_",joys[i],"_1_red.dput"))
    schmertout2_red[[i]] <- dget(paste0("output/schmert_",joys[i],"_2_red.dput"))
  }
}

#ASFRs
schmertadat <- jccadat %>% left_join(alldata2 %>% filter(par==0) %>% select(ctry,Cohort,Age,x)) %>%
  mutate(pred=NA,sd=NA)
for (i in 1:length(joys)) {
  for (j in ctry_joy_df %>% filter(joy==joys[i]) %>% pull(ctry)) {
    j2 <- j
    ind <- which(schmertadat$ctry==j & schmertadat$joy==joys[i])
    if (ctry_joy_df %>% filter(joy==joys[i],ctry==j) %>% pull(red) == F) {
      tmp <- schmertout1[[i]][[j2]]
    } 
    else { tmp <- schmertout1_red[[i]][[j2]] }
    schmertadat$pred[ind] <- tmp$mu %>% as.data.frame() %>%
      mutate(Age=agerange) %>%
      pivot_longer(-Age,names_to="Cohort",values_to="pred") %>%
      mutate(Cohort=as.numeric(Cohort)) %>%
      mutate(joy=joys[i],ctry=j) %>% arrange(Cohort,Age) %>%
      select(pred) %>% unlist()
    schmertadat$sd[ind] <- tmp$sd %>% as.data.frame() %>%
      mutate(Age=agerange) %>%
      pivot_longer(-Age,names_to="Cohort",values_to="sd") %>%
      mutate(Cohort=as.numeric(Cohort)) %>%
      mutate(joy=joys[i],ctry=j) %>% arrange(Cohort,Age) %>%
      select(sd) %>% unlist()
  }
}
schmertadat %<>% mutate(
  lwr95=as.vector(qnorm(0.025,pred,sd)),
  upr95=as.vector(qnorm(0.975,pred,sd)),
  lwr90=as.vector(qnorm(0.05,pred,sd)),
  upr90=as.vector(qnorm(0.95,pred,sd)),
  lwr50=as.vector(qnorm(0.25,pred,sd)),
  upr50=as.vector(qnorm(0.75,pred,sd)))

#Calculate CRPS
schmertadat$crps <- NA
ind <- which(schmertadat$Year>schmertadat$joy & !is.na(schmertadat$x))
schmertadat$crps[ind] <- CRPS_normal_func(schmertadat$x[ind],schmertadat$pred[ind],schmertadat$sd[ind])

#Calculate LogS
schmertadat$logs <- NA
schmertadat$logs[ind] <- -log(dnorm(schmertadat$x[ind],schmertadat$pred[ind],schmertadat$sd[ind])) 

#CFRs
#Version 1: normal
schmertcdat <- jccdat %>% left_join(obscfr2) %>% mutate(pred=NA,sd=NA)
for (i in 1:length(joys)) {
  for (j in ctry_joy_df %>% filter(joy==joys[i]) %>% pull(ctry)) {
    j2 <- j
    ind <- which(schmertcdat$ctry==j & schmertcdat$joy==joys[i])
    if (ctry_joy_df %>% filter(joy==joys[i],ctry==j) %>% pull(red) == F) {
      tmp <- schmertout1[[i]][[j2]]
    } 
    else { tmp <- schmertout1_red[[i]][[j2]] }
    schmertcdat$pred[ind] <- tmp$mu.cfr
    schmertcdat$sd[ind] <- tmp$sd.cfr
  }
}
schmertcdat %<>% mutate(
  lwr95=as.vector(qnorm(0.025,pred,sd)),
  upr95=as.vector(qnorm(0.975,pred,sd)),
  lwr90=as.vector(qnorm(0.05,pred,sd)),
  upr90=as.vector(qnorm(0.95,pred,sd)),
  lwr50=as.vector(qnorm(0.25,pred,sd)),
  upr50=as.vector(qnorm(0.75,pred,sd)))

#Version 2: treating observed ASFRs as having no error (only used for quantitative comparisons with other methods)
schmertjcdat <- jccdat %>% left_join(obscfr) %>% left_join(obscfr2) %>% mutate(pred=NA,sd=NA)
for (i in 1:length(joys)) {
  for (j in ctry_joy_df %>% filter(joy==joys[i]) %>% pull(ctry)) {
    j2 <- j
    ind <- which(schmertjcdat$ctry==j & schmertjcdat$joy==joys[i])
    if (ctry_joy_df %>% filter(joy==joys[i],ctry==j) %>% pull(red) == F) {
      tmp <- schmertout2[[i]][[j2]]
    } 
    else { tmp <- schmertout2_red[[i]][[j2]] }
    schmertjcdat$pred[ind] <- tmp$mu.cfr+schmertjcdat$jx[ind]
    schmertjcdat$sd[ind] <- tmp$sd.cfr
  }
}

#Join with uncertainty intervals from version 1
schmertjcdat %<>% left_join(schmertcdat %>% select(Cohort,ctry,joy,starts_with("lwr"),starts_with("upr")))

#Calculate CRPS
schmertjcdat$crps <- NA
ind <- which(schmertjcdat$Cohort+max(agerange)>schmertjcdat$joy & !is.na(schmertjcdat$x))
schmertjcdat$crps[ind] <- CRPS_normal_func(schmertjcdat$x[ind],schmertjcdat$pred[ind],schmertjcdat$sd[ind])

#Calculate LogS
schmertjcdat$logs <- NA
schmertjcdat$logs[ind] <- -log(dnorm(schmertjcdat$x[ind],schmertjcdat$pred[ind],schmertjcdat$sd[ind]))


##UN model (model of Ševčíková et al.)
#Fit UN model
source("r/5a_fit_UN_model.r")

#Prepare age-period exposures
source("r/5a_prepare_age_period_exp.r")

#Compute ASFR and CFR samples and scoring rules
unjcdat <- unadat <- data.frame() #to store results
for (joy2 in joys) { #for each JOY
  print(joy2)
  #Read in ASFR samples
  file <- paste0("output/bayesTFR_asfr_sample_all_",joy2,".txt")
  tmp <- read_tsv(file)
  #Compute CFRs observed up to the year JOY
  tmp2p <- jccadat %>% filter(joy == joy2, Year<=joy2) %>% left_join(alldata2 %>% filter(par==0)) %>%
    group_by(ctry,Cohort) %>% mutate(n=n()) %>% filter(n<30) %>% summarise(jx=sum(x))
  #Convert ASFR samples into long data frame (only for the year JOY+1 onwards)
  tmp %<>% mutate(ctry=countrycode(sourcevar=country_code,origin="iso3n",destination="iso3c")) %>%
    filter(ctry %in% unique(tmp2p$ctry)) %>%
    rename(Age=age,Year=year,Cohort=coh) %>%
    select(Age, ctry, Year, Cohort, paste(1:500)) %>%
    pivot_longer(-c("Age","Year","Cohort","ctry"),names_to="iter",values_to="ASFR") %>%
    mutate(iter=as.numeric(iter))
  #Join with age-period exposures
  tmp %<>% left_join(APcomb %>% rename(ctry=Code)) %>% mutate(b=N*ASFR)
  #Add on a duplicate of each row for the previous cohort (for each age-period combination,
  #we observe the cohort a+y and a+y-1)
  tmp %<>% rbind(tmp %>% mutate(Cohort=Cohort-1))
  #Halve birth counts so total births preserved, and retain only age-cohort combinations with both observations
  tmp %<>% select(-c("ASFR","N")) %>% mutate(b=0.5*b) %>%
    group_by(ctry,Age,Cohort,iter) %>% mutate(n=n()) %>% filter(n==2) %>%
    filter(Cohort %in% unique(tmp2p$Cohort))
  #Aggregate births by country, age and cohort
  tmp %<>% group_by(ctry,Age,Cohort,iter) %>% summarise(b=sum(b))
  #Construct year column
  tmp %<>% mutate(Year=Age+Cohort)
  #Join with age-cohort exposures, compute age-cohort ASFRs
  tmp %<>% left_join(ACexpmat2) %>% mutate(ASFR=b/N) %>% select(-c("b","N"))
  #Retain ASFR samples for later use
  tmpa <- tmp
  
  #Construct temporary ASFR summary data frame
  tmp_asum <- tmpa %>% group_by(ctry,Cohort,Age) %>%
    summarise(med=median(ASFR),
              mean=mean(ASFR),
              sd=sd(ASFR),
              lwr95=quantile(ASFR,p=0.025),
              upr95=quantile(ASFR,p=0.975),
              lwr90=quantile(ASFR,p=0.05),
              upr90=quantile(ASFR,p=0.95),
              lwr50=quantile(ASFR,p=0.25),
              upr50=quantile(ASFR,p=0.75))
  
  #Construct temporary ASFR results data frame to extract
  unadat_tmp <- jccadat %>% filter(joy==joy2) %>%
    left_join(alldata2 %>% filter(par==0) %>% select(ctry,Cohort,Age,x)) %>%
    left_join(tmp_asum %>% rename(pred=med))
  
  #Compute sum of ASFRs
  tmp %<>% group_by(ctry,Cohort,iter) %>% summarise(ASFR=sum(ASFR))
  #Compute overall CFR prediction
  tmp %<>% left_join(tmp2p) %>% mutate(CFR=ASFR+jx) %>% select(ctry,Cohort,iter,CFR)
  
  #Construct temporary CFR summary data frame
  tmp_sum <- tmp %>% group_by(ctry,Cohort) %>%
    summarise(med=median(CFR),
              mean=mean(CFR),
              sd=sd(CFR),
              lwr95=quantile(CFR,p=0.025),
              upr95=quantile(CFR,p=0.975),
              lwr90=quantile(CFR,p=0.05),
              upr90=quantile(CFR,p=0.95),
              lwr50=quantile(CFR,p=0.25),
              upr50=quantile(CFR,p=0.75))
  
  #Construct temporary CFR results data frame to extract
  unjcdat_tmp <- jccdat %>% filter(joy==joy2) %>%
    left_join(obscfr) %>% left_join(obscfr2) %>%
    left_join(tmp_sum %>% rename(pred=med))
  
  #CRPS for CFR
  unjcdat_tmp$crps <- NA
  ind <- which(!is.na(unjcdat_tmp$pred) & !is.na(unjcdat_tmp$x))
  Ns <- 500
  for (i in ind) {
    cfr_sam <- tmp %>% filter(ctry==unjcdat_tmp$ctry[i],Cohort==unjcdat_tmp$Cohort[i]) %>% pull(CFR)
    unjcdat_tmp$crps[i] <- (2/(Ns^2))*sum((cfr_sam - unjcdat_tmp$x[i])*(Ns*as.numeric(rep(unjcdat_tmp$x[i], Ns) < cfr_sam) - (1:Ns) + 0.5))
  }
  
  #LogS for CFR
  unjcdat_tmp$logs <- NA
  unjcdat_tmp$logs[ind] <- -log(dnorm(unjcdat_tmp$x[ind],unjcdat_tmp$mean[ind],unjcdat_tmp$sd[ind]))
  
  #CRPS for ASFR
  unadat_tmp$crps <- NA
  ind <- which(!is.na(unadat_tmp$pred) & !is.na(unadat_tmp$x))
  Ns <- 500
  for (i in ind) {
    asfr_sam <- tmpa %>% filter(ctry==unadat_tmp$ctry[i],Cohort==unadat_tmp$Cohort[i],Age==unadat_tmp$Age[i]) %>% pull(ASFR)
    unadat_tmp$crps[i] <- (2/(Ns^2))*sum((asfr_sam - unadat_tmp$x[i])*(Ns*as.numeric(rep(unadat_tmp$x[i], Ns) < asfr_sam) - (1:Ns) + 0.5))
  }
  
  #LogS for ASFR
  unadat_tmp$logs <- NA
  unadat_tmp$logs[ind] <- -log(dnorm(unadat_tmp$x[ind],unadat_tmp$mean[ind],unadat_tmp$sd[ind]))
  
  #Extract results
  unjcdat <- rbind(unjcdat,unjcdat_tmp)
  unadat <- rbind(unadat,unadat_tmp)
}
save(unjcdat,unadat,file="output/UN_sumstat.RData")
#load("output/UN_sumstat.RData")


##Proposed model (aggregate and parity-specific)
ctryvec <- rep(ctry_joy_df$ctry,2)
ctrylvec <- rep(ctrys_long[ctryvec])
mods <- c("Aggregate","Parity-specific")
modvec <- rep(rep(mods),each=nrow(ctry_joy_df))
joyvec <- rep(ctry_joy_df$joy,2)
redvec <- rep(ctry_joy_df$red,2)
files <- paste0("mod_",substr(tolower(modvec),1,1),"_",ctryvec,"_",joyvec)
for (i in 1:nrow(ctry_joy_df)) {
  if (ctry_joy_df$red[i]==T) files[c(i,i+nrow(ctry_joy_df))] <- paste0(files[c(i,i+nrow(ctry_joy_df))],"_red")
}
files <- paste0(files,".RData")

genunc <- gen_unc2(files,agerange,joyvec,modvec,ctryvec,redvec)
save(genunc,file="output/prop_sumstat.RData")
#load("output/prop_sumstat.RData")

#Process output
#ASFRs
propadat <- genunc$asfr
propadat$pred <- propadat$med
propadat %<>% left_join(data.frame(mod=factor(1:length(files)),
                                   Model=modvec,ctry=ctryvec,joy=joyvec),multiple="all")
propadat %<>% select(colnames(propadat)[c(2,3,19,20,4,14,17,7:13,15:16,18)]) %>%
  rename("Age"=age,"Cohort"=coh,"Year"=yr)

#CFRs
propcdat <- genunc$cfr
propcdat %<>% left_join(data.frame(mod=factor(1:length(files)),
                                    Model=modvec,ctry=ctryvec,joy=joyvec),multiple="all")
propjcdat <- genunc$jcfr
propjcdat$pred <- propjcdat$med2
propjcdat %<>% left_join(data.frame(mod=factor(1:length(files)),
                                    Model=modvec,ctry=ctryvec,joy=joyvec),multiple="all") %>%
  select(-starts_with("lwr"),-starts_with("upr")) %>%
  left_join(propcdat %>% select(mod,coh,starts_with("lwr"),starts_with("upr")))
propjcdat %<>% select(colnames(propjcdat)[c(14,2,15,12,7,5,16:21,10:11,13)]) %>%
  rename("Cohort"=coh,"x"=ocfr)
