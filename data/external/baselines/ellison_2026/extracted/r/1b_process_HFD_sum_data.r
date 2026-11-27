###Processing HFD data - summary fertility measures
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
HFDdat %<>% filter(min_age==min(agerange),max_age==max(agerange),n_age==Na)


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

#Remove country-cohort combinations with unreliable estimates
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

#Compute scaled standard deviations for comparison
HFDdat5$sd/c(mean(HFDdat3$value[HFDdat3$name=="PPR01"],na.rm=T),
             mean(HFDdat3$value[HFDdat3$name=="PPR12"],na.rm=T),
             mean(HFDdat3$value[HFDdat3$name=="PPR23"],na.rm=T),
             mean(HFDdat3$value[HFDdat3$name=="PPR3p4p"],na.rm=T))
HFDdat6$sd/mean(HFDdat3$value[HFDdat3$name=="CFR"],na.rm=T)


##Plot original and differenced PPRs
supp.ppr <- c("CFR","PPR 0 -> 1","PPR 1 -> 2","PPR 2 -> 3","PPR 3+ -> 4+")
names(supp.ppr) <- c("CFR",paste0("PPR",c("01","12","23","3p4p")))

png(file=paste0("plots/HFD_sum_original_",joy,"_red.png"),width=27.5,height=10,units="cm",res=400,family="Calibri")
print(HFDdat3 %>% filter(substr(name,1,3)=="PPR" | name=="CFR") %>%
  ggplot(aes(x=Cohort,y=value,col=Code)) +
  facet_wrap(~name,scales="free_y",labeller = labeller(name=supp.ppr),nrow=1) +
  labs(y="Value",col="Country") +
  geom_line(show.legend = F) + theme_bw() +
  coord_cartesian(xlim=c(min_coh-4,max_coh+4)) +
  guides(color=guide_legend(nrow=1,byrow=T)) +
  scale_x_continuous(breaks=seq(1905,2005,10),minor_breaks=setdiff(1900:2000,seq(1905,2005,10))) +
  geom_dl(aes(label = Code), method = list(fontfamily="Calibri",cex=0.75,dl.combine("first.points", "last.points"))) +
  theme(legend.position = "bottom",text=element_text("Calibri")))
dev.off()

png(file=paste0("plots/HFD_sum_original_2row_",joy,"_red.png"),width=15,height=15,units="cm",res=400,family="Calibri")
print(HFDdat3 %>% filter(substr(name,1,3)=="PPR" | name=="CFR") %>%
        ggplot(aes(x=Cohort,y=value,col=Code)) +
        facet_wrap(~name,scales="free_y",nrow=2,labeller = labeller(name=supp.ppr)) +
        labs(y="Value",col="Country") +
        geom_line(show.legend = F) + theme_bw() +
        coord_cartesian(xlim=c(min_coh-5,max_coh+5)) +
        guides(color=guide_legend(nrow=1,byrow=T)) +
        scale_x_continuous(breaks=seq(1905,2005,10),minor_breaks=setdiff(1900:2000,seq(1905,2005,10))) +
        geom_dl(aes(label = Code), method = list(fontfamily="Calibri",cex=0.75,dl.combine("first.points", "last.points"))) +
        theme(legend.position = "bottom",text=element_text("Calibri")))
dev.off()

png(file=paste0("plots/HFD_sum_differenced_",joy,"_red.png"),width=27.5,height=10,units="cm",res=400,family="Calibri")
print(HFDdat4 %>% filter(substr(name,1,3)=="PPR" | name == "CFR") %>%
  ggplot(aes(x=Cohort,y=value,col=Code)) +
  facet_wrap(~name,labeller = labeller(name=supp.ppr),nrow=1) +
  labs(y="Value",col="Country") +
  geom_line(show.legend = F) + theme_bw() +
  coord_cartesian(xlim=c(min_coh-4,max_coh+4),ylim=c(-0.1,0.1)) +
  guides(color=guide_legend(nrow=1,byrow=T)) +
  scale_x_continuous(breaks=seq(1905,2005,10),minor_breaks=setdiff(1900:2000,seq(1905,2005,10))) +
  geom_dl(aes(label = Code), method = list(fontfamily="Calibri",cex=0.75,dl.combine("first.points", "last.points"))) +
  geom_hline(yintercept=0) +
  theme(legend.position = "bottom",text=element_text("Calibri")))
dev.off()
