###Processing HFD data
##Extract HFD data
#Extract "asfrVHbo.txt" from asfr.zip
filenames <- as.character(unzip("data/asfr.zip", list = TRUE)$Name)
unzip("data/asfr.zip", filenames[which(filenames=="Files/zip_w/asfrVHbo.txt")], exdir="data", junkpaths=T)


##Process data
#Read in HFD data
HFDdat <- read.table("data/asfrVHbo.txt",skip=2,as.is=T, header=T, na.strings = ".") %>% as_tibble()
#Change Age column to numeric
HFDdat$Age <- as.numeric(HFDdat$Age)
#Remove rows with missing ASFRs
HFDdat <- HFDdat[which(!is.na(HFDdat$ASFR)),]
#Only include ages within the agerange considered
HFDdat %<>% filter(Age %in% agerange)
#For each country and cohort, find minimum and maximum ages included, and the length of the agerange
HFDdat %<>% group_by(Code,Cohort) %>% mutate(min_age=min(Age),max_age=max(Age),n_age=length(Age)) %>% ungroup()
#Only select complete cohorts (when age_min = start of agerange, age_max = end of agerange,
#and n_age = length of agerange)
HFDdat %<>% filter(min_age==min(agerange),max_age==max(agerange),n_age==length(agerange))


##Prepare uniform distributions
resSMunif <- list()
for (j in 1:length(joys)) {
  resSMunif[[j]] <- list()
  joy <- joys[j]
  ##Compute cohort summary fertility measures
  HFDdat2 <- HFDdat %>% group_by(Code,Cohort) %>%
    summarise(CFR=sum(ASFR),
              CFR1=sum(ASFR1),
              CFR2=sum(ASFR2),
              CFR3=sum(ASFR3),
              CFR4=sum(ASFR4),
              CFR5=sum(ASFR5p),
              CFR3p=sum(ASFR3+ASFR4+ASFR5p),
              CFR4p=sum(ASFR4+ASFR5p),
              PPR01=CFR1,
              PPR12=CFR2/CFR1,
              PPR23=CFR3/CFR2,
              PPR3p4p=CFR4p/CFR3p)
  #Select cohort range for HFD analysis
  min_coh <- 1945
  max_coh <- joy-max(agerange)
  HFDdat2 %<>% filter(Cohort >= min_coh, Cohort <= max_coh)
  #Select countries that have at least 2 complete cohorts (to provide at least one cohort difference)
  code_inc <- HFDdat2 %>% group_by(Code) %>% summarise(x=length(unique(Cohort))) %>% filter(x>=2) %>%
    ungroup() %>% select(Code) %>% as.data.frame()
  HFDdat2 %<>% filter(Code %in% code_inc$Code)
  
  #Convert into long data frame
  HFDdat3 <- HFDdat2 %>% pivot_longer(cols=starts_with(c("CFR","PPR")))
  
  #Remove countries with unreliable estimates
  for (i in names(coh_iss)) {
    HFDdat3$value[which(HFDdat3$Code==i & HFDdat3$Cohort %in% coh_iss[[i]])] <- NA
  }
  
  #Subset to countries considered for fitting
  HFDdat3 %<>% filter(Code %in% ctrys)
  
  #Compute differenced values
  HFDdat4 <- HFDdat3 %>% group_by(Code,name) %>%
    reframe(Cohort=Cohort, value=value-lag(value)) %>% filter(!is.na(value)) %>% ungroup()
  
  ##Compute standard deviation of differenced PPRs
  HFDdat5 <- HFDdat4 %>% filter(substr(name,1,3)=="PPR") %>% group_by(name) %>%
    summarise(sd=sd(value))
  HFDdat6 <- HFDdat4 %>% filter(name=="CFR") %>% group_by(name) %>%
    summarise(sd=sd(value))
  
  ##Compute country-specific SDs
  supp.ppr <- c("CFR","PPR 0 -> 1","PPR 1 -> 2","PPR 2 -> 3","PPR 3+ -> 4+")
  names(supp.ppr) <- c("CFR",paste0("PPR",c("01","12","23","3p4p")))
  
  #Must have at least 5 differenced observations per measure to be included (otherwise
  #we are not fitting to the country anyway)
  HFDunif1 <- HFDdat4 %>% filter(substr(name,1,3)=="PPR" | name=="CFR") %>%
    filter(!is.na(value)) %>%
    group_by(Code,name) %>% mutate(n=n()) %>%
    filter(n>=5) %>%
    group_by(Code,name) %>% summarise(sd=sd(value)) %>%
    left_join(rbind(HFDdat5,HFDdat6) %>% rename(sdt=sd)) %>%
    group_by(name) %>% arrange(name,desc(sd)) %>% mutate(ind=1:n())  
  
  #Pick countries to fit to
  ctrys_fit_joy <- lapply(joy_sel_lst,function(x) joy %in% x) %>% unlist()
  resSMunif[[j]]$ctry <- HFDunif1 %>% filter(name=="CFR") %>% ungroup() %>%
    filter(Code %in% ctrys[which(ctrys_fit_joy==T)]) %>%
    mutate(min=min(sd),max=max(sd),diff=sd-sdt, mindiff=min(abs(diff))) %>%
    filter(min==sd | max==sd | mindiff==abs(diff)) %>% arrange(sd)
  
  #Assign lower and upper bounds for uniform priors
  resSMunif[[j]]$bounds <- HFDunif1 %>% group_by(name) %>% 
    mutate(min=min(sd),max=max(sd)) %>%
    filter(min==sd | max==sd) %>% arrange(name,sd)
  
  png(file=paste0("plots/HFD_sum_sd_ctry_",joy,"_SMunif.png"),width=15,height=15,units="cm",res=400,family="Calibri")
  print({ HFDunif1 %>%
      ggplot(aes(x=ind,y=sd,color=Code)) +
      geom_hline(aes(yintercept=sdt),linetype=2) +
      geom_point() +
      facet_wrap(~name,labeller=labeller(name=supp.ppr)) +
      coord_cartesian(ylim=c(0,0.035)) +
      scale_x_continuous(breaks=0:100,minor_breaks=NULL) +
      geom_text(aes(ind,sd+0.003,label=Code),
                angle=90, size = 2, vjust=0.4, show.legend = F) +
      labs(x="",y="SD",color="Country") +
      theme_bw() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) })
  dev.off()
}

##Save output
save(resSMunif,file="output/resSMunif.RData")
