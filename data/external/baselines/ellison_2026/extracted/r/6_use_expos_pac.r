###Assessing effect of exposure approximation on forecasts

##Read in register-based exposure estimates
expos  <-  read.table("data/exposRRpac.txt",skip=2,as.is=T, header=T) %>% as_tibble()


##Extract available years for each country
expos %<>% filter(Code %in% ctrys_fit)
ctry_yr_pac <- list()
ctrys_pac <- numeric()
for (i in ctrys_fit) { #for each country
  yrs <- unique(expos$Year[expos$Code==i]) #Extract available years 
  ctry_yr_pac[[i]] <- yrs #assign to list
  #Check if we have at least 36 years to fit at least the reduced cohort model
  if (length(yrs)>=36 & all(diff(yrs)==1)) ctrys_pac <- c(ctrys_pac,i)
}
ctry_yr_pac[ctrys_pac]

#Years needed for each JOY - 40 cohorts
yfunc <- function(joy,red) {
  if (red==F) cohrange <- (joy-max(agerange)-10):(joy-min(agerange))
  if (red==T) cohrange <- (joy-max(agerange)-5):(joy-min(agerange))
  yrs <- c(cohrange+15,max(cohrange)+15+1)
  c(min(yrs),max(yrs))
}
yrs_full <- sapply(joys,yfunc,F)
yrs_red <- sapply(joys,yfunc,T)

#Identify available JOYs for each country
joy_sel_lst_pac <- red_sel_lst_pac <- list()
for (i in ctrys_pac) {
  joy_sel <- joy_sel_lst[[i]]
  red_sel <- red_sel_lst[[i]]
  for (j in 1:length(joy_sel)) {
    yrs <- yfunc(joy_sel[j],red_sel[j])
    if (min(ctry_yr_pac[[i]])<=yrs[1] & max(ctry_yr_pac[[i]])>=yrs[2]) {
      joy_sel_lst_pac[[i]] <- c(joy_sel_lst_pac[[i]],joy_sel[j])
      red_sel_lst_pac[[i]] <- c(red_sel_lst_pac[[i]],red_sel[j])
    }
    else
      if (red_sel[j] == F) {
      yrs <- yfunc(joy_sel[j],T)
      if (min(ctry_yr_pac[[i]])<=yrs[1] & max(ctry_yr_pac[[i]])>=yrs[2]) {
        joy_sel_lst_pac[[i]] <- c(joy_sel_lst_pac[[i]],joy_sel[j])
        red_sel_lst_pac[[i]] <- c(red_sel_lst_pac[[i]],T)
      }
    }
  }
}


##Compute rates
for (ctry in ctrys_pac) {
  print(ctry)
  joy_sel <- joy_sel_lst_pac[[ctry]]
  red_sel <- red_sel_lst_pac[[ctry]]
  for (i in 1:length(joy_sel)) {
    #Set JOY
    joy <- joy_sel[i]; print(joy)
    red <- red_sel[i]
    #Set cohort range (if possible, 40 cohorts are considered: 11 complete cohorts and age 15 only 
    #observed in the last cohort; otherwise, 35 cohorts are considered, with only 6 complete cohorts)
    if (red==F) cohrange <- (joy-max(agerange)-10):(joy-min(agerange))
    if (red==T) cohrange <- (joy-max(agerange)-5):(joy-min(agerange))
    
    #Compare rate estimates and save rates
    source("r/6_prepare_data_pac_plot.r")
  } 
}


##Prepare data for fitting models, fit models, plot model output
for (ctry in ctrys_pac) {
  print(ctry)
  joy_sel <- joy_sel_lst_pac[[ctry]]
  red_sel <- red_sel_lst_pac[[ctry]]
  for (i in 1:length(joy_sel)) {
    #Set JOY
    joy <- joy_sel[i]; print(joy)
    red <- red_sel[i]
    #Set cohort range (if possible, 40 cohorts are considered: 11 complete cohorts and age 15 only 
    #observed in the last cohort; otherwise, 35 cohorts are considered, with only 6 complete cohorts)
    if (red==F) cohrange <- (joy-max(agerange)-10):(joy-min(agerange))
    if (red==T) cohrange <- (joy-max(agerange)-5):(joy-min(agerange))
    
    #Fit original parity-specific model if it doesn't already exist
    if (!(red==F & file.exists(paste0("output/mod_p_",ctry,"_",joy,".RData")) |
        red==T & file.exists(paste0("output/mod_p_",ctry,"_",joy,"_red.RData")))) {
      source("r/2_prepare_stan_data.r")
      source("r/3_fit_mod_p.r")
    }
    
    #Prepare stan data
    if (!(red==F & file.exists(paste0("data/standata_",ctry,"_",joy,"_pac.RData")) |
          red==T & file.exists(paste0("data/standata_",ctry,"_",joy,"_pac_red.RData"))))
    source("r/6_prepare_stan_data_pac.r")
    
    #Model 2
    if ((red==F & file.exists(paste0("output/mod_p_pac2_",ctry,"_",joy,".RData")) |
          red==T & file.exists(paste0("output/mod_p_pac2_",ctry,"_",joy,"_red.RData")))) {
      source("r/6_fit_mod_p_pac2.r")
      source("r/6_plot_model_output_ctry_joy_pac2.r")
    }
    
    #Model 3
    if ((red==F & file.exists(paste0("output/mod_p_pac3_",ctry,"_",joy,".RData")) |
          red==T & file.exists(paste0("output/mod_p_pac3_",ctry,"_",joy,"_red.RData")))) {
      source("r/6_fit_mod_p_pac3.r")
      source("r/6_plot_model_output_ctry_joy_pac3.r")
    }
  } 
}


##Summarise model output
comp_time <- rhat <- list()
for (ctry in ctrys_pac) {
  comp_time[[ctry]] <- rhat[[ctry]] <- list()
  print(ctry)
  joy_sel <- joy_sel_lst_pac[[ctry]]
  red_sel <- red_sel_lst_pac[[ctry]]
  for (i in 1:length(joy_sel)) {
    comp_time[[ctry]][[i]] <- rhat[[ctry]][[i]] <- list()
    #Set JOY
    joy <- joy_sel[i]; print(joy)
    red <- red_sel[i]
    #Set cohort range (if possible, 40 cohorts are considered: 11 complete cohorts and age 15 only 
    #observed in the last cohort; otherwise, 35 cohorts are considered, with only 6 complete cohorts)
    if (red==F) cohrange <- (joy-max(agerange)-10):(joy-min(agerange))
    if (red==T) cohrange <- (joy-max(agerange)-5):(joy-min(agerange))
    
    #Model 2
    if (red==F) file <- paste0("output/mod_p_pac2_",ctry,"_",joy,".RData")
    if (red==T) file <- paste0("output/mod_p_pac2_",ctry,"_",joy,"_red.RData")
    if (file.exists(file)) {
      load(file)
      comp_time[[ctry]][[i]][[2]] <- get_elapsed_time(stanout) %>% sum()/60
      a <- summary(stanout)$summary
      rhat[[ctry]][[i]][[2]] <- range(a[,"Rhat"],na.rm=T)
    }
    
    #Model 3
    if (red==F) file <- paste0("output/mod_p_pac3_",ctry,"_",joy,".RData")
    if (red==T) file <- paste0("output/mod_p_pac3_",ctry,"_",joy,"_red.RData")
    if (file.exists(file)) {
      load(file)
      comp_time[[ctry]][[i]][[3]] <- get_elapsed_time(stanout) %>% sum()/60
      a <- summary(stanout)$summary
      rhat[[ctry]][[i]][[3]] <- range(a[,"Rhat"],na.rm=T)
    }

  }
}

comp_time
rhat


##Compare forecast performance
#Combine datasets
alldata2 <- data.frame()
for (i in ctrys_pac) {
  load(paste0("data/",i,"_rates.RData"))
  alldata2 <- rbind(alldata2,data.frame(ctry=i,alldata))
}
alldata2 %<>% mutate(Year=Cohort+Age) #add year column

#Setup
nmod <- 3
sumdat <- data.frame(ctry=rep(rep(ctrys_pac,lapply(joy_sel_lst_pac,length)),each=nmod),
                     joy=rep(unlist(joy_sel_lst_pac),each=3),
                     red=rep(unlist(red_sel_lst_pac),each=3))

ctryvec <- sumdat$ctry
ctrylvec <- rep(ctrys_long[ctryvec])
mods <- paste0(rep(c("Parity-specific, "),nrow(sumdat)),
               rep(c("Original","Variant 1","Variant 2"),nrow(sumdat)/nmod))
modvec <- mods
joyvec <- sumdat$joy
redvec <- sumdat$red
regexp <- rep(c(F,T,T),nrow(sumdat)/nmod)

#Prepare file names of fitted objects
modt <- rep(c("Original","Variant 1","Variant 2"),nrow(sumdat)/nmod)
files <- numeric()
for (i in 1:length(modvec)) {
  files[i] <- paste0("mod_",substr(tolower(modvec[i]),1,1))
  if (modt[i]=="Variant 1") files[i] <- paste0(files[i],"_pac2")
  if (modt[i]=="Variant 2") files[i] <- paste0(files[i],"_pac3")
  files[i] <- paste0(files[i],"_",ctryvec[i],"_",joyvec[i])
  if (redvec[i]==T) files[i] <- paste0(files[i],"_red")
}
files <- paste0(files,".RData")

genunc <- gen_unc2(files,agerange,joyvec,modvec,ctryvec,redvec)
save(genunc,file="output/prop_sumstat_pac.RData")
#load("output/prop_sumstat_pac.RData")


##CFRs
#Extract summary statistics
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

#Prepare for comparison
combcdat <- propjcdat %>% 
  mutate(Model=factor(Model,levels=unique(modvec))) %>%
  arrange(Model,joy,ctry,Cohort)
combcdat2 <- combcdat %>% filter(Cohort+max(agerange)>joy,!is.na(x))
combcdat2$err <- combcdat2$x-combcdat2$pred

#Create tables of overall summary statistics
cfunc <- function() {
  combcdat2 %>%
    group_by(Model) %>%  
    summarise(
      MAE=mean(abs(err)),
      MAPE=100*mean(abs((err)/x)),
      RMSE=sqrt(mean(err^2)),
      Bias=mean(-err),
      Rel_bias=100*mean(-err/x),
      cov90=100*mean(x>=lwr90 & x<=upr90),
      cov50=100*mean(x>=lwr50 & x<=upr50),
      CRPS=mean(crps)) %>%
    flextable() %>%
    set_formatter(
      MAE = function(x) formatC(x, digits = 3, format = 'f'),
      MAPE = function(x) formatC(x, digits = 1, format = 'f'),
      RMSE = function(x) formatC(x, digits = 3, format = 'f'),
      Bias = function(x) formatC(x, digits = 3, format = 'f'),
      Rel_bias = function(x) formatC(x, digits = 1, format = 'f'),
      cov90 = function(x) formatC(x, digits = 0, format = 'f'),
      cov50 = function(x) formatC(x, digits = 0, format = 'f'),
      CRPS = function(x) formatC(x, digits = 3, format = 'f')) %>%
    bold(~ min(MAPE) == MAPE, ~ MAPE) %>%
    bold(~ min(RMSE) == RMSE, ~ RMSE) %>%
    bold(~ min(MAE) == MAE, ~ MAE) %>%
    bold(~ min(abs(Bias)) == abs(Bias), ~ Bias) %>%
    bold(~ min(abs(Rel_bias)) == abs(Rel_bias), ~ Rel_bias) %>%
    bold(~ min(abs(cov90-90),na.rm=T) == abs(cov90-90), ~ cov90) %>%
    bold(~ min(abs(cov50-50),na.rm=T) == abs(cov50-50), ~ cov50) %>%
    bold(~ min(CRPS) == CRPS, ~ CRPS) %>%
    font(part = 'all', fontname = 'Calibri') %>%
    fontsize(part = 'all', size = 10) %>%
    bold(part = 'header') %>%
    align(part = 'header', align = 'center') %>%
    set_caption(paste0("CFR summary statistics by model: Parity-specific, all countries, all JOYs"))
}

tmp2 <- cfunc()

#Setup for table
tfiles <- c('tables/CFR_pac.docx')
sect_properties <- prop_section(
  page_size = page_size(
    orient = "portrait",
    width = 8.3, height = 11.7
  ),
  type = "continuous",
  page_margins = page_mar()
)

#Save tables to file
save_as_docx(tmp2,pr_section=sect_properties, path = tfiles)

#Prepare dataset to plot CRPS
combacdat2 <- data.frame(type="CFR",
                         combcdat2 %>% select(Model,ctry,Cohort,joy,crps) %>%
                           group_by(Model,ctry,joy) %>% summarise(crps=mean(crps)))

combacdat2 %<>% left_join(data.frame(Model=unique(modvec),
                                     Model2=c("Original","Variant 1","Variant 2")))
combacdat2 %<>% left_join(data.frame(ctry=ctrys,ctryl=ctrys_long))

#Combine with coverage
combacdat3 <- data.frame(type="CFR",
                         combcdat2 %>%
                           group_by(Model,ctry,joy) %>%
                           summarise(crps=mean(crps),
                                     cov90=mean(x>=lwr90 & x<=upr90),
                                     cov50=mean(x>=lwr50 & x<=upr50)))

combacdat3 %<>% left_join(data.frame(Model=unique(modvec),
                                     Model3=c("Original","Variant 1","Variant 2")))
combacdat3 %<>% left_join(data.frame(ctry=ctrys,ctryl=ctrys_long))

combacdat3 %<>% pivot_longer(c("crps","cov90","cov50"))

#Prepare colours for plotting
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#Save plot of CRPS, 50% and 90% coverage - version 1
png(file=paste0("plots/CRPS_cov50_90_ctry_pac_CFR.png"),width=10,height=15,units="cm",res=400,family="Calibri")
supp.labv <- c("Average CRPS","Coverage of 50% PI","Coverage of 90% PI"); names(supp.labv) <- c("crps","cov50","cov90")
print({combacdat3 %>%
    ggplot(aes(x=ctry,y=value)) +
    facet_grid(name~paste0("JOY = ",joy),scales="free_y",labeller=labeller(name=supp.labv)) +
    labs(y="Value",x="Country",color="Model") +
    scale_color_manual(values=c(1,gg_color_hue(3))) +
    geom_hline(data=data.frame(name=c("cov50","cov90"),yintercept=c(0.5,0.9)),
               aes(yintercept=yintercept),linetype=2) +
    geom_errorbar(data=combacdat3 %>% filter(Model3=="Original"),
                  aes(ymin=value,ymax=value),width=0.5,color=1) +
    geom_point(aes(color=Model3)) +
    theme_bw() +
    theme(legend.position="bottom")})
dev.off()

#Save plot of CRPS, 50% and 90% coverage - version 2
png(file=paste0("plots/CRPS_cov50_90_ctry_pac_CFR_2.png"),width=10,height=15,units="cm",res=400,family="Calibri")
supp.labv <- c("Average CRPS","Coverage of 50% PI","Coverage of 90% PI"); names(supp.labv) <- c("crps","cov50","cov90")
print({combacdat3 %>%
    ggplot(aes(x=factor(joy),y=value)) +
    facet_grid(name~ctryl,scales="free_y",labeller=labeller(name=supp.labv)) +
    labs(y="Value",x="JOY",color="Model") +
    scale_color_manual(values=c(1,gg_color_hue(3))) +
    geom_hline(data=data.frame(name=c("cov50","cov90"),yintercept=c(0.5,0.9)),
               aes(yintercept=yintercept),linetype=2) +
    geom_errorbar(data=combacdat3 %>% filter(Model3=="Original"),
                  aes(ymin=value,ymax=value),width=0.5,color=1) +
    geom_point(aes(color=Model3)) +
    theme_bw() +
    theme(legend.position="bottom")})
dev.off()


##ASFRs
#Extract summary statistics
propadat <- genunc$asfr
propadat$pred <- propadat$med
propadat %<>% left_join(data.frame(mod=factor(1:length(files)),
                                   Model=modvec,ctry=ctryvec,joy=joyvec),multiple="all")
propadat %<>% select(colnames(propadat)[c(2,3,19,20,4,14,17,7:13,15:16,18)]) %>%
  rename("Age"=age,"Cohort"=coh,"Year"=yr)

#Prepare for comparison
combadat <- propadat %>% 
  mutate(Model=factor(Model,levels=unique(modvec))) %>%
  arrange(Model,joy,ctry,Cohort,Age)
combadat2 <- combadat %>% filter(Year>joy,!is.na(x))
combadat2$err <- combadat2$x-combadat2$pred

#Create tables of overall summary statistics
afunc <- function() {
  tmp2 <- combadat2 %>%
    group_by(Model) %>%  
    summarise(
      MAE=mean(abs(err)),
      MAPE=100*mean(abs((err)/x)),
      RMSE=sqrt(mean(err^2)),
      Bias=mean(-err),
      Rel_bias=100*mean(-err/x),
      cov90=100*mean(x>=lwr90 & x<=upr90),
      cov50=100*mean(x>=lwr50 & x<=upr50),
      CRPS=mean(crps)) %>%
    flextable() %>%
    set_formatter(
      MAE = function(x) formatC(x, digits = 3, format = 'f'),
      MAPE = function(x) formatC(x, digits = 1, format = 'f'),
      RMSE = function(x) formatC(x, digits = 3, format = 'f'),
      Bias = function(x) formatC(x, digits = 3, format = 'f'),
      Rel_bias = function(x) formatC(x, digits = 1, format = 'f'),
      cov90 = function(x) formatC(x, digits = 0, format = 'f'),
      cov50 = function(x) formatC(x, digits = 0, format = 'f'),
      CRPS = function(x) formatC(x, digits = 3, format = 'f')) %>%
    bold(~ min(MAPE) == MAPE, ~ MAPE) %>%
    bold(~ min(RMSE) == RMSE, ~ RMSE) %>%
    bold(~ min(MAE) == MAE, ~ MAE) %>%
    bold(~ min(abs(Bias)) == abs(Bias), ~ Bias) %>%
    bold(~ min(abs(Rel_bias)) == abs(Rel_bias), ~ Rel_bias) %>%
    bold(~ min(abs(cov90-90),na.rm=T) == abs(cov90-90), ~ cov90) %>%
    bold(~ min(abs(cov50-50),na.rm=T) == abs(cov50-50), ~ cov50) %>%
    bold(~ min(CRPS) == CRPS, ~ CRPS) %>%
    font(part = 'all', fontname = 'Calibri') %>%
    fontsize(part = 'all', size = 10) %>%
    bold(part = 'header') %>%
    align(part = 'header', align = 'center') %>%
    set_caption(paste0("ASFR summary statistics by model: Parity-specific, all countries, all JOYs"))
}

tmp2 <- afunc()

#Setup for table
tfiles <- c('tables/ASFR_pac.docx')
sect_properties <- prop_section(
  page_size = page_size(
    orient = "portrait",
    width = 8.3, height = 11.7
  ),
  type = "continuous",
  page_margins = page_mar()
)

#Save tables to file
save_as_docx(tmp2, pr_section=sect_properties, path = tfiles)

#Prepare dataset to plot CRPS
combacdat2 <- data.frame(type="ASFR",
                         combadat2 %>% select(Model,ctry,Cohort,joy,crps) %>%
                           group_by(Model,ctry,joy) %>% summarise(crps=mean(crps)))

combacdat2 %<>% left_join(data.frame(Model=unique(modvec),
                                     Model3=c("Original","Variant 1","Variant 2")))
combacdat2 %<>% left_join(data.frame(ctry=ctrys,ctryl=ctrys_long))

#Combine with coverage
combacdat3 <- data.frame(type="ASFR",
                         combadat2 %>%
                           group_by(Model,ctry,joy) %>%
                           summarise(crps=mean(crps),
                                     cov90=mean(x>=lwr90 & x<=upr90),
                                     cov50=mean(x>=lwr50 & x<=upr50)))

combacdat3 %<>% left_join(data.frame(Model=unique(modvec),
                                     Model3=c("Original","Variant 1","Variant 2")))
combacdat3 %<>% left_join(data.frame(ctry=ctrys,ctryl=ctrys_long))

combacdat3 %<>% pivot_longer(c("crps","cov90","cov50"))

#Save plot of CRPS, 50% and 90% coverage - version 1
png(file=paste0("plots/CRPS_cov50_90_ctry_pac_ASFR.png"),width=10,height=15,units="cm",res=400,family="Calibri")
supp.labv <- c("Average CRPS","Coverage of 50% PI","Coverage of 90% PI"); names(supp.labv) <- c("crps","cov50","cov90")
print({combacdat3 %>%
    ggplot(aes(x=ctry,y=value)) +
    facet_grid(name~paste0("JOY = ",joy),scales="free_y",labeller=labeller(name=supp.labv)) +
    labs(y="Value",x="Country",color="Model") +
    scale_color_manual(values=c(1,gg_color_hue(3))) +
    geom_hline(data=data.frame(name=c("cov50","cov90"),yintercept=c(0.5,0.9)),
               aes(yintercept=yintercept),linetype=2) +
    geom_errorbar(data=combacdat3 %>% filter(Model3=="Original"),
                  aes(ymin=value,ymax=value),width=0.5,color=1) +
    geom_point(aes(color=Model3)) +
    theme_bw() +
    theme(legend.position="bottom")})
dev.off()

#Save plot of CRPS, 50% and 90% coverage - version 2
png(file=paste0("plots/CRPS_cov50_90_ctry_pac_ASFR_2.png"),width=10,height=15,units="cm",res=400,family="Calibri")
supp.labv <- c("Average CRPS","Coverage of 50% PI","Coverage of 90% PI"); names(supp.labv) <- c("crps","cov50","cov90")
print({combacdat3 %>%
    ggplot(aes(x=factor(joy),y=value)) +
    facet_grid(name~ctryl,scales="free_y",labeller=labeller(name=supp.labv)) +
    labs(y="Value",x="JOY",color="Model") +
    scale_color_manual(values=c(1,gg_color_hue(3))) +
    geom_hline(data=data.frame(name=c("cov50","cov90"),yintercept=c(0.5,0.9)),
               aes(yintercept=yintercept),linetype=2) +
    geom_errorbar(data=combacdat3 %>% filter(Model3=="Original"),
                  aes(ymin=value,ymax=value),width=0.5,color=1) +
    geom_point(aes(color=Model3)) +
    theme_bw() +
    theme(legend.position="bottom")})
dev.off()
