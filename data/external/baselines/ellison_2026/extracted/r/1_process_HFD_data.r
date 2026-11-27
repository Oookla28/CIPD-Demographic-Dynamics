###Process HFD/HMD data
##Extract HFD data
#Extract births data ('birthsVH.txt' and 'birthsVHBo.txt' from 'births.zip')
filenames <- as.character(unzip("data/births.zip", list = TRUE)$Name)
unzip("data/births.zip", filenames[which(filenames=="Files/zip_w/birthsVH.txt")], exdir="data", junkpaths=T)
unzip("data/births.zip", filenames[which(filenames=="Files/zip_w/birthsVHbo.txt")], exdir="data", junkpaths=T)

#Extract exposure data ('exposVH.txt' from 'exposure.zip')
filenames <- as.character(unzip("data/exposure.zip", list = TRUE)$Name)
unzip("data/exposure.zip", filenames[which(filenames=="Files/zip_w/exposVH.txt")], exdir="data", junkpaths=T)


##Process aggregate data
#Read in births data from HFD
birth  <-  read.table("data/birthsVH.txt",skip=2,as.is=T, header=T) %>% as_tibble()

#Read in exposure data from HFD
expos  <-  read.table("data/exposVH.txt",skip=2,as.is=T, header=T) %>% as_tibble()

#Restrict to country and reproductive age range
birth %<>% filter(Code==ctry, Age %in% 15:44) %>%
  mutate(Age=as.numeric(Age), Total=as.numeric(Total))
expos %<>% filter(Code==ctry, Age %in% 15:44) %>%
  mutate(Exposure=as.numeric(Exposure))

#Round births and exposures
birth$Total <- round(birth$Total)
expos$Exposure <- round(expos$Exposure)

#Combine births and exposure data
dat <- birth %>% left_join(expos) %>% rename(b=Total,N=Exposure)
dat$ASFR <- dat$b/dat$N

#Remove rows with missing births/exposures
dat %<>% filter(!is.na(ASFR))
  
#Assign dat to datagg
datagg <- dat


##Process order-specific data
#Read in births by birth order data from HFD
birth  <-  read.table("data/birthsVHbo.txt",skip=2, as.is=T, header=T) %>% as_tibble()

#Restrict to country and reproductive age range
birth %<>% filter(Code==ctry,Age %in% 15:44)

#Make age and birth columns rounded numeric variables
for (col in c("Age","Total",paste0("B",c(1:4,"5p"))))
  birth[[col]] <- round(as.numeric(birth[[col]]))

#Recalculate birth totals
birth$Total <- rowSums(birth[,paste0("B",c(1:4,"5p"))])

#Combine births and exposure data
dat <- birth %>% left_join(expos) %>% rename(N=Exposure)


##Compute parity-specific rate estimates
#Setup
coh_sel <- dat %>% filter(Age==15, !is.na(Total)) %>% pull(Cohort)
dat %<>% filter(Cohort %in% coh_sel)
dat$E5p <- dat$E4 <- dat$E3 <- dat$E2 <- dat$E1 <- NA
dat$r5p <- dat$r4 <- dat$r3 <- dat$r2 <- dat$r1 <- NA
dat$Esum <- NA

#Loop through cohorts
for (i in unique(dat$Cohort)) {
  tmp <- dat[dat$Cohort==i,]
  ind <- match(paste0("E",c(1:4,"5p")),colnames(tmp))
  tmp[1,ind] <- c(tmp[1,"N"],rep(0,4))
  if (nrow(tmp)>1) {
  for (j in 2:nrow(tmp)) {
    tmp[j,ind[1]] <- tmp[j,"N"]*(tmp[j-1,"E1"]-tmp[j-1,"B1"])/tmp[j-1,"N"]
    tmp[j,ind[2]] <- tmp[j,"N"]*(tmp[j-1,"E2"]+tmp[j-1,"B1"]-tmp[j-1,"B2"])/tmp[j-1,"N"]
    tmp[j,ind[3]] <- tmp[j,"N"]*(tmp[j-1,"E3"]+tmp[j-1,"B2"]-tmp[j-1,"B3"])/tmp[j-1,"N"]
    tmp[j,ind[4]] <- tmp[j,"N"]*(tmp[j-1,"E4"]+tmp[j-1,"B3"]-tmp[j-1,"B4"])/tmp[j-1,"N"]
    tmp[j,ind[5]] <- tmp[j,"N"]*(tmp[j-1,"E5p"]+tmp[j-1,"B4"])/tmp[j-1,"N"]
    tmp[j,ind[1:4]] <- round(tmp[j,ind[1:4]])
    tmp[j,ind[5]] <- tmp[j,"N"] - sum(tmp[j,ind[1:4]])
  }
  }
  tmp$Esum <- apply(tmp[,ind],1,sum)
  tmp$r1 <- tmp$B1/tmp$E1
  tmp$r2 <- tmp$B2/tmp$E2
  tmp$r3 <- tmp$B3/tmp$E3
  tmp$r4 <- tmp$B4/tmp$E4
  tmp$r5p <- tmp$B5p/tmp$E5p
  dat[dat$Cohort==i,] <- tmp
}

#Extract results
datt <- rbind(data.frame(par=0,datagg[,2:3],B=datagg$b,E=datagg$N),
              data.frame(par=1,dat[,2:3],B=dat$B1,E=dat$E1),
              data.frame(par=2,dat[,2:3],B=dat$B2,E=dat$E2),
              data.frame(par=3,dat[,2:3],B=dat$B3,E=dat$E3),
              data.frame(par=4,dat[,2:3],B=dat$B4,E=dat$E4),
              data.frame(par=5,dat[,2:3],B=dat$B5p,E=dat$E5p))


##Process WPP exposures for future age-cohort combinations (after last year of HFD data)
#Read in WPP estimates and projections
load("data/WPP_pop.RData")
#Subset to country
ACall %<>% filter(Country==ctry)
#Subset to years after last observed year of HFD data
ACproj <- ACall %>% mutate(Year=age+Cohort) %>%
  filter(Year > max(datagg$Age + datagg$Cohort))
#Rename age column
ACproj %<>% rename(Age=age)
#Round population counts
ACproj %<>% mutate(N=round(N))


##Process and save data
#Combine parities 3 and above
dattb <- datt
dattb$par[dattb$par>4] <- 4
dattb %<>% group_by(par,Cohort,Age) %>% summarise(N=sum(E),n=sum(B),x=n/N) %>% ungroup()

#Remove invalid rates
dattb <- dattb[which(dattb$N>0 & dattb$n>=0 & dattb$n <= dattb$N),]
alldata <- dattb

#Combine observed and projected exposures
#Extract midyear population estimates from 'alldata'
alldata2 <- alldata %>% filter(par==0)
#Combine with WPP exposures
allexp <- rbind(alldata2[,2:4],ACproj[,c("Cohort","Age","N")])
#Create data frame with age and cohort (Cohort) columns
ACexpmat <- data.frame(expand.grid(Age=15:44,Cohort=unique(allexp$Cohort)))
#Add year column
ACexpmat$Year <- ACexpmat$Age + ACexpmat$Cohort
#Add population estimates (N) column
ACexpmat %<>% left_join(allexp)


##Save data
save(alldata,ACexpmat,file=paste0("data/",ctry,"_rates.RData"))
