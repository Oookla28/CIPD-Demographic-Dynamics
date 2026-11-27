###Prepare data for modelling
##Read in rates
load(paste0("data/",ctry,"_rates.RData"))


##Modelling setup
Na <- length(agerange) # number of ages (30)
Nc <- length(cohrange) # number of cohorts (40 or 35)
#Create data frame of all age-cohort combinations
dat <- expand.grid(Age=agerange,Cohort=cohrange)
Nf <- nrow(dat) # number of age-cohort combinations (30 x 40 or 35 = 1200 or 1050)
#Extract age-cohort exposures for full surface
ACexp <- ACexpmat$N[ACexpmat$Cohort %in% cohrange & ACexpmat$Age %in% agerange]

#Number of age-cohort combinations for each parity (Nf for ASFRs/parity 0, then
#Nc fewer for each successive parity, as we start from age 16 for parity 1, age 17
#for parity 2, and age 18 for parity 3+)
Nf1 <- Nf
Nf2 <- Nf-Nc
Nf3 <- Nf-2*Nc
Nf4 <- Nf-3*Nc

#Create indicators for age-cohort combinations omitted and included for parities 1, 2 and 3+
indA2 <- which(dat$Age<16)
indA2b <- setdiff(1:Nf,indA2)
indA3 <- which(dat$Age<17)
indA3b <- setdiff(1:Nf,indA3)
indA4 <- which(dat$Age<18)
indA4b <- setdiff(1:Nf,indA4)


##Prepare data
alldata$Year <- alldata$Age + alldata$Cohort #add year column
#All (ages 15-44)
data0 <- alldata[alldata$Age %in% agerange & alldata$Cohort %in% cohrange & alldata$par==0,]
indp0 <- which(!is.na(left_join(data.frame(expand.grid(Age=agerange,Cohort=cohrange)),data0)$x))
N0 <- nrow(data0)
succ0 <- data0$n
tot0 <- data0$N

#Parity 0 (ages 15-44)
data1 <- alldata[alldata$Age %in% agerange & alldata$Cohort %in% cohrange & alldata$par==1,]
indp1 <- which(!is.na(left_join(data.frame(expand.grid(Age=agerange,Cohort=cohrange)),data1)$par))
N1 <- nrow(data1)
succ1 <- data1$n

#Parity 1 (ages 16-44)
data2 <- alldata[alldata$Age %in% agerange[-c(1)] & alldata$Cohort %in% cohrange & alldata$par==2,]
indp2 <- which(!is.na(left_join(data.frame(expand.grid(Age=agerange,Cohort=cohrange)),data2)$par))
N2 <- nrow(data2)
succ2 <- data2$n

#Parity 2 (ages 17-44)
data3 <- alldata[alldata$Age %in% agerange[-c(1,2)] & alldata$Cohort %in% cohrange & alldata$par==3,]
indp3 <- which(!is.na(left_join(data.frame(expand.grid(Age=agerange,Cohort=cohrange)),data3)$par))
N3 <- nrow(data3)
succ3 <- data3$n

#Parity 3+ (ages 18-44)
data4 <- alldata[alldata$Age %in% agerange[-c(1,2,3)] & alldata$Cohort %in% cohrange & alldata$par==4,]
indp4 <- which(!is.na(left_join(data.frame(expand.grid(Age=agerange,Cohort=cohrange)),data4)$par))
N4 <- nrow(data4)
succ4 <- data4$n


##Additional data preparation to incorporate observed exposures
expos  <-  read.table("data/exposRRpac.txt",skip=2,as.is=T, header=T) %>% as_tibble()
#Restrict to country and reproductive age range
expos %<>% filter(Code==ctry, Age %in% agerange)
#Make into a long data frame
expos %<>% pivot_longer(-c("Code","Year","Age"),names_to="par",values_to="N")
#Make par a numeric variable
expos %<>% mutate(par=as.numeric(substr(par,2,2))+1)
#Combine parities 3 and above
expos %<>% mutate(par=ifelse(par>4,4,par))
expos %<>% group_by(par,Year,Age) %>% summarise(N=sum(N))
#Add cohort column (y-a, the younger cohort)
expos %<>% mutate(Cohort=Year-Age)
#Duplicate expos, with the duplication having the older cohort (y-(a+1))
expos <- rbind(expos,expos %>% mutate(Cohort=Cohort-1))
#Halve the exposures
expos %<>% mutate(N=0.5*N)
#Add a column to count the number of copies of each age-cohort combination
expos %<>% group_by(par,Age,Cohort) %>% mutate(ind=n())
#Remove age-cohort combinations that are only partially observed
expos %<>% filter(ind==2)
#Sum exposures across country, age and cohort
expos %<>% group_by(par,Age,Cohort) %>% summarise(N=sum(N))
#Arrange by cohort
expos %<>% arrange(par,Cohort,Age)
#Filter to included cohorts
expos %<>% filter(Cohort %in% cohrange)
#Merge in observed age-cohort total exposures
expos %<>% left_join(ACexpmat %>% rename(Ntot=N))
#Add column with implied total from adjusted exposures
expos %<>% group_by(Age,Cohort) %>% mutate(Ntot2=sum(N))
#Compare AC totals
#plot(expos$Ntot,expos$Ntot2);abline(0,1)
#Rescale exposures to match those in ACexpmat
expos %<>% mutate(N=N*Ntot/Ntot2)
#Check totals match
expos %<>% mutate(Ntot2b=sum(N))
#Compare AC totals
#plot(expos$Ntot,expos$Ntot2b);abline(0,1)
#Return to a wide data frame
expos %<>% pivot_wider(id_cols=c("Age","Cohort"),names_from=par,values_from=N)
#Create objects for model fitting
tot_pac <- expos[,-c(1:2)] %>% as.matrix()
tot1 <- data1 %>% select(-N) %>% left_join(expos %>% select(Age,Cohort,`1`))
tot2 <- data2 %>% select(-N) %>% left_join(expos %>% select(Age,Cohort,`2`))
tot3 <- data3 %>% select(-N) %>% left_join(expos %>% select(Age,Cohort,`3`))
tot4 <- data4 %>% select(-N) %>% left_join(expos %>% select(Age,Cohort,`4`))

##Prepare rates
rdat <- rbind(data0 %>% mutate(N2=N) %>% select(-N) %>% rename(N=N2) %>% mutate(x2=x),
              tot1 %>% rename(N=`1`) %>% mutate(x2=n/N),
              tot2 %>% rename(N=`2`) %>% mutate(x2=n/N),
              tot3 %>% rename(N=`3`) %>% mutate(x2=n/N),
              tot4 %>% rename(N=`4`) %>% mutate(x2=n/N))

#Save rates for later use
rdat2 <- rdat %>% select(-x) %>% rename(x=x2)
filename <- ifelse(red==F,paste0("data/",ctry,"_rates_pac_",joy,".RData"),
                   paste0("data/",ctry,"_rates_pac_",joy,"_red.RData"))
save(rdat2,file=filename)

##Plot rates
supp.labs1 <- c("ASFR","Parity 0","Parity 1","Parity 2","Parity 3+")
names(supp.labs1) <- c(0,1,2,3,4)
supp.labs2 <- c("Cumulation-based exposures","Register-based exposures")
names(supp.labs2) <- c("x","x2")

plotname <- ifelse(red==F,paste0("plots/",ctry,"_rates_pac_",joy,".png"),
                   paste0("plots/",ctry,"_rates_pac_",joy,"_red.png"))
png(file=plotname,width=20,height=10,units="cm",res=400,family="Calibri")
print({rdat %>% filter(par>0) %>% select(-c("N","n","Year")) %>%
    pivot_longer(-c("par","Cohort","Age")) %>%
  ggplot(aes(x=Age,y=value)) +
  geom_line(aes(group=Cohort,color=Cohort),lwd=0.3) +
  facet_grid(name~par) +
  scale_color_gradientn(colours=viridis(100)[100:1], guide = guide_colorbar(barheight = 15,frame.colour="black",ticks.colour="black"), breaks = seq(1800,2000,5), limits=range(cohrange)) +
  scale_y_continuous(breaks=seq(0,0.5,0.1)) +
  scale_x_continuous(breaks=seq(15,45,5), minor_breaks=setdiff(15:45,seq(15,45,5)), expand=c(0.02,0.02)) +
  coord_cartesian(ylim=c(0,0.5)) +
  theme_bw() + theme(text = element_text("Calibri")) +
  facet_grid(name~par,labeller=labeller(par=supp.labs1,name=supp.labs2))})
dev.off()
