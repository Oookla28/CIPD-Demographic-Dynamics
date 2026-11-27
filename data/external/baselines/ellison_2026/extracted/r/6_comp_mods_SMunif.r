###Compare models with uniform prior for summary measures
load("output/resSMunif.RData")
ctrys_bf <- lapply(resSMunif, function(x) x$ctry$Code) %>% unlist()


##Prepare colours for plotting
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


##Setup
#Combine datasets
alldata2 <- data.frame()
for (i in ctrys_fit) {
  load(paste0("data/",i,"_rates.RData"))
  alldata2 <- rbind(alldata2,data.frame(ctry=i,alldata))
}
alldata2 %<>% mutate(Year=Cohort+Age) #add year column

#Construct minimum cohort data frame
ctry_joy_df <- data.frame(
  ctry=substr(joy_sel_lst %>% unlist() %>% names(),1,3),
  joy=joy_sel_lst %>% unlist(),
  red=red_sel_lst %>% unlist())
ctry_joy_df %<>% mutate(mincoh=ifelse(red==T,joy-max(agerange)-5,joy-max(agerange)-10))
ctry_joy_df2 <- data.frame()
for (i in 1:length(joys)) {
  ctry_joy_df2 <- rbind(ctry_joy_df2,
                        ctry_joy_df %>% filter(joy==joys[i],
                                               ctry %in% resSMunif[[i]]$ctry$Code))
}


##Proposed model (aggregate and parity-specific)
#Setup
ctryvec <- rep(ctrys_bf,2*2)
ctrylvec <- rep(ctrys_long[ctryvec])
mods <- paste0(rep(c("Aggregate, ","Parity-specific, "),each=2*length(ctrys_bf)),
               rep(rep(c("Original (fixed SD)",
                     "Uniform prior on SD"),each=length(ctrys_bf)),2))
modvec <- mods
joyvec <- rep(rep(joys,each=3),2*2)
redvec <- rep(ctry_joy_df2$red[match(paste0(ctrys_bf,rep(joys,each=3)),
                                     apply(ctry_joy_df2[,1:2],1,function(x) paste0(x[1],x[2])))],2*2)

#Prepare file names of fitted objects
modt <- rep(rep(c("Original","Uniform"),each=length(ctrys_bf)),2)
files <- numeric()
for (i in 1:length(modvec)) {
  files[i] <- paste0("mod_",substr(tolower(modvec[i]),1,1))
  if (modt[i]=="Uniform") files[i] <- paste0(files[i],"_SMunif")
  files[i] <- paste0(files[i],"_",ctryvec[i],"_",joyvec[i])
  if (redvec[i]==T) files[i] <- paste0(files[i],"_red")
}
files <- paste0(files,".RData")

genunc <- gen_unc2(files,agerange,joyvec,modvec,ctryvec,redvec)
save(genunc,file="output/prop_sumstat_SMunif.RData")
#load("output/prop_sumstat_SMunif.RData")


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
mods2 <- c("Aggregate","Parity-specific")
cfunc <- function(i) {
  combcdat2 %>% filter(substr(Model,1,1)==substr(mods2,1,1)[i]) %>%
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
  set_caption(paste0("CFR summary statistics by model: ",mods2[i],", all countries, all JOYs"))
}

tmp2 <- cfunc(1)
tmp3 <- cfunc(2)

#Setup for table
tfiles <- c('tables/CFR_SMunif.docx')
sect_properties <- prop_section(
  page_size = page_size(
    orient = "portrait",
    width = 8.3, height = 11.7
  ),
  type = "continuous",
  page_margins = page_mar()
)

#Save tables to file
save_as_docx(tmp2,tmp3,pr_section=sect_properties, path = tfiles)

#Prepare dataset to plot CRPS
combacdat2 <- data.frame(type="CFR",
                         combcdat2 %>% select(Model,ctry,Cohort,joy,crps) %>%
                           group_by(Model,ctry) %>% summarise(crps=mean(crps)))

combacdat2 %<>% left_join(data.frame(Model=unique(modvec),
                                     Model2=rep(c("Aggregate","Parity-specific"),each=2),
                                     Model3=rep(c("Original","Uniform"),2)))
combacdat2 %<>% left_join(data.frame(ctry=ctrys,ctryl=ctrys_long))

#Combine with coverage
combacdat3 <- data.frame(type="CFR",
                         combcdat2 %>%
                           group_by(Model,ctry) %>%
                           summarise(crps=mean(crps),
                                     cov90=mean(x>=lwr90 & x<=upr90),
                                     cov50=mean(x>=lwr50 & x<=upr50)))

combacdat3 %<>% left_join(data.frame(Model=unique(modvec),
                                     Model2=rep(c("Aggregate","Parity-specific"),each=2),
                                     Model3=rep(c("Original","Uniform prior"),2)))
combacdat3 %<>% left_join(data.frame(ctry=ctrys,ctryl=ctrys_long))

combacdat3 %<>% pivot_longer(c("crps","cov90","cov50"))

#Save plot of CRPS, 50% and 90% coverage
png(file=paste0("plots/CRPS_cov50_90_ctry_SMunif_CFR.png"),width=12.5,height=17.5,units="cm",res=400,family="Calibri")
supp.labv <- c("Average CRPS","Coverage of 50% PI","Coverage of 90% PI"); names(supp.labv) <- c("crps","cov50","cov90")
print({combacdat3 %>%
    ggplot(aes(x=ctry,y=value)) +
    facet_grid(name~Model2,scales="free_y",labeller=labeller(name=supp.labv)) +
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
mods2 <- c("Aggregate","Parity-specific")
afunc <- function(i) {
  tmp2 <- combadat2 %>% filter(substr(Model,1,1)==substr(mods2[i],1,1)) %>%
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
    set_caption(paste0("ASFR summary statistics by model: ",mods2[i],", all countries, all JOYs"))
}

tmp2 <- afunc(1)
tmp3 <- afunc(2)

#Setup for table
tfiles <- c('tables/ASFR_SMunif.docx')
sect_properties <- prop_section(
  page_size = page_size(
    orient = "portrait",
    width = 8.3, height = 11.7
  ),
  type = "continuous",
  page_margins = page_mar()
)

#Save tables to file
save_as_docx(tmp2, tmp3, pr_section=sect_properties, path = tfiles)

#Prepare dataset to plot CRPS
combacdat2 <- data.frame(type="ASFR",
                         combadat2 %>% select(Model,ctry,Cohort,joy,crps) %>%
                           group_by(Model,ctry) %>% summarise(crps=mean(crps)))

combacdat2 %<>% left_join(data.frame(Model=unique(modvec),
                                     Model2=rep(c("Aggregate","Parity-specific"),each=2),
                                     Model3=rep(c("Original","Uniform"),2)))
combacdat2 %<>% left_join(data.frame(ctry=ctrys,ctryl=ctrys_long))

#Combine with coverage
combacdat3 <- data.frame(type="ASFR",
                         combadat2 %>%
                           group_by(Model,ctry) %>%
                           summarise(crps=mean(crps),
                                     cov90=mean(x>=lwr90 & x<=upr90),
                                     cov50=mean(x>=lwr50 & x<=upr50)))

combacdat3 %<>% left_join(data.frame(Model=unique(modvec),
                                     Model2=rep(c("Aggregate","Parity-specific"),each=2),
                                     Model3=rep(c("Original","Uniform prior"),2)))
combacdat3 %<>% left_join(data.frame(ctry=ctrys,ctryl=ctrys_long))

combacdat3 %<>% pivot_longer(c("crps","cov90","cov50"))

#Save plot of CRPS, 50% and 90% coverage
png(file=paste0("plots/CRPS_cov50_90_ctry_SMunif_ASFR.png"),width=12.5,height=17.5,units="cm",res=400,family="Calibri")
supp.labv <- c("Average CRPS","Coverage of 50% PI","Coverage of 90% PI"); names(supp.labv) <- c("crps","cov50","cov90")
print({combacdat3 %>%
    ggplot(aes(x=ctry,y=value)) +
    facet_grid(name~Model2,scales="free_y",labeller=labeller(name=supp.labv)) +
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
