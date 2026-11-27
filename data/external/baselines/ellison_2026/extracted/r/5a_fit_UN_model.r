###UN fertility projections
##Extract HFD data
#Extract "asfrRR.txt" from asfr.zip
filenames <- as.character(unzip("data/asfr.zip", list = TRUE)$Name)
unzip("data/asfr.zip", filenames[which(filenames=="Files/zip_w/asfrRR.txt")], exdir="data", junkpaths=T)
#Extract "exposRR.txt" from expos.zip
filenames <- as.character(unzip("data/exposure.zip", list = TRUE)$Name)
unzip("data/exposure.zip", filenames[which(filenames=="Files/zip_w/exposRR.txt")], exdir="data", junkpaths=T)


##Read in HFD data
HFDdat <- read.table("data/asfrRR.txt",skip=2,as.is=T, header=T, na.strings = ".")
#Restrict to agerange
HFDdat %<>% filter(Age %in% agerange)
#Check no missing data
sum(is.na(HFDdat)) == 0
#Check all 30 ages observed in each country-year combination
HFDdat %>% group_by(Code,Year) %>% mutate(n=n()) %>% pull(n) %>% table()


##Compute TFR
tdat <- HFDdat %>% group_by(Code,Year) %>% summarise(TFR=sum(ASFR))
#Edit country codes
tdat$Code[tdat$Code=="FRATNP"] <- "FRA"
tdat$Code[tdat$Code=="DEUTNP"] <- "DEU"
tdat$Code[tdat$Code=="GBR_NP"] <- "GBR"
#Add country names and ISO codes
tdat %<>% left_join(codelist %>% select(country.name.en,genc3c,iso3n) %>%
                      rename(name=country.name.en,Code=genc3c,country_code=iso3n))
tdat %>% filter(is.na(name)) %>% pull(Code) %>% unique() #countries with missing codes
tdat %<>% filter(!is.na(name)) #remove countries that do not have a code
#Check all countries included in UNlocations dataset
data(UNlocations)
setdiff(unique(tdat$country_code),UNlocations$country_code)

#Compute PASFRs
HFDdat$Code[HFDdat$Code=="FRATNP"] <- "FRA"
HFDdat$Code[HFDdat$Code=="DEUTNP"] <- "DEU"
HFDdat$Code[HFDdat$Code=="GBR_NP"] <- "GBR"
pdat <- HFDdat %>% left_join(tdat) %>% filter(!is.na(name))
pdat %<>% group_by(Code,Year) %>% mutate(PASFR=100*ASFR/TFR)
#Add in zero columns for ages 10-14, 45-54
pdat <- rbind(
  pdat,
  pdat %>% filter(Age %in% paste(15:29)) %>% mutate(ASFR=0,PASFR=0) %>%
                group_by(Code,Year) %>% mutate(Age=paste(c(10:14,45:54)))) %>%
  arrange(Code,Year,Age)


##Fit TFR and ASFR models
for (joy in joys) { #for each JOY
  #Restrict data to the years (joy-30):joy, and at least the two most recent TFR values observed
  tdat2 <- tdat %>% group_by(Code) %>%
    mutate(miny=min(Year),maxy=max(Year)) %>%
    filter(maxy>=joy) %>%
    filter(miny<=joy-1) %>%
    filter(Year>=(joy-30),Year<=joy) %>%
    arrange(Year,Code)
  #Check all countries in ctrys_fit included in tdat2
  setdiff(ctrys_fit,unique(tdat2$Code))
  ctrycodes <- unique(tdat2$country_code)
  ctrynames <- unique(tdat2$Code)
  
  #Rearrange and save dataset
  sample_tfr <- tdat2 %>%
    pivot_wider(id_cols=c("country_code","name"),
                values_from="TFR",
                names_from="Year")
  write_tsv(sample_tfr,file=paste0("data/sample_tfr_",joy,".txt"))

  #Run phase 2 model
  m2 <- run.tfr.mcmc(nr.chains=1,
                     iter=60000,
                     output.dir=paste0("bayesTFR.output.",joy),
                     thin=1,
                     annual=TRUE,
                     ar.phase2=TRUE,
                     uncertainty=FALSE,
                     start.year=joy-30,
                     present.year=joy,
                     use.wpp.data=FALSE,
                     my.tfr.file = paste0("data/sample_tfr_",joy,".txt"),
                     replace.output = TRUE)
  
  #Inspect output
  tfr.dir <- paste0("bayesTFR.output.",joy)
  m2 <- get.tfr.mcmc(tfr.dir)
  m2$mcmc.list
  tfr.partraces.plot(m2,sim.dir=tfr.dir)
 
  #Check there are the correct number of samples
  tmpdir <- paste0(tfr.dir,"/mc1")
  filenames <- list.files(tmpdir)
  filenames <- setdiff(filenames,"bayesTFR.mcmc.rda")
  nrowdat <- data.frame()
  for (i in 1:length(filenames)) {
    tmp <- read.table(paste0(tmpdir,"/",filenames[i]), quote="\"", comment.char="")
    nrowdat <- rbind(nrowdat,data.frame(filenames=filenames[i],nrow=nrow(tmp)))
  }
  nrowdat_sel <- nrowdat %>% filter(nrow!=60000)
  nrowdat_sel #should be empty

  #Run phase 3 model
  m3 <- run.tfr3.mcmc(sim.dir = tfr.dir,
                      nr.chains = 1,
                      iter = 60000, 
                      thin = 1,
                      replace.output = TRUE)
  
  #Inspect output
  m3 <- get.tfr3.mcmc(tfr.dir)
  m3$mcmc.list
  tfr3.partraces.plot(m3)
  
  #Extract phase 3 countries
  ph3ctrycode <- m3$meta$regions$country_code[m3$meta$id_phase3]
  ph3ctryname <- m3$meta$regions$country_name[m3$meta$id_phase3]
  
  #Generate TFR predictions
  tfr.pred <- tfr.predict(
    sim.dir = tfr.dir,
    start.year = joy+1,
    end.year = joy+30, 
    burnin = 30000,
    use.tfr3 = T,
    burnin3 = 30000,
    nr.traj=500,
    use.correlation = T,
    replace.output=T)
  
  #Inspect output
  tfr.pred <- get.tfr.prediction(tfr.dir)
  
  #Extract and save TFR trajectories for later use
  convert.tfr.trajectories(tfr.dir,"all")
  
  #Prepare PASFRs
  pdat2 <- pdat %>%
    filter(Code %in% tdat2$Code) %>%
    filter(Year %in% (joy-30):joy) %>%
    pivot_wider(id_cols=c("country_code","name","Age"),
                values_from="PASFR",
                names_from="Year") %>% rename(age=Age)
  write_tsv(pdat2,file=paste0("data/obs_pasfr_",joy,".txt"))

  #Prepare pattern
  patdat <- tdat2 %>% filter(Year==joy) %>% ungroup() %>%
    select(country_code,name,Code) %>%
    rename(country=name) %>%
    mutate(PasfrNorm="GlobalNorm")
  #Implement weighting procedure as in Bohk-Ewald et al. 2018
  wdat <- pdat %>% filter(Year==joy,Code %in% tdat2$Code) %>%
    group_by(country_code,Code) %>%
    summarise(TFR=mean(TFR),
              MAB=sum((as.numeric(Age)+0.5)*ASFR)/sum(ASFR)) %>%
    ungroup()
  wdat %<>% mutate(minTFR=min(TFR),maxMAB=max(MAB),
                   dTFR=(max(TFR)-TFR)/(max(TFR)-min(TFR)),
                   dMAB=(MAB-min(MAB))/(max(MAB)-min(MAB)))
  wdat %<>% left_join(data.frame(country_code=ph3ctrycode,Ip3=1)) %>%
    mutate(Ip3=ifelse(is.na(Ip3),0,Ip3))
  wdat %<>% mutate(wt=0.2*dTFR+0.7*dMAB+0.1*Ip3)
  wdat %<>% arrange(desc(wt))
  n20p <- floor(nrow(wdat)*0.2)
  wdat2 <- wdat %>% head(n20p)
  
  #Produce ASFR projections
  dir.create(paste0(tfr.dir,"/pasfr_traj"))
  for (ctry in ctrycodes) {
    patdat$PasfrGlobalNorm <- 0
    wdat_tmp <- wdat2 %>% filter(wt>=wdat$wt[which(wdat$country_code==ctry)])
    if (nrow(wdat_tmp)==1) wdat_tmp <- wdat[1:2,]
    patdat$PasfrGlobalNorm[match(wdat_tmp$country_code,patdat$country_code)] <- 1
    write_tsv(patdat,file=paste0("data/pattern_",ctry,"_",joy,".txt"))  
  sample_pasfr <- 
    project.pasfr.traj(inputs=list(tfr.file=paste0(tfr.dir,"/predictions/ascii_trajectories.csv"),
                                   pasfr=paste0("data/obs_pasfr_",joy,".txt"),
                                   patterns=paste0("data/pattern_",ctry,"_",joy,".txt")),
                       countries=ctry,
                       present.year=joy,
                       end.year=joy+30,
                       out.file.name=paste0(tfr.dir,"/pasfr_traj/pasfr_",ctry,".txt"),
                       annual=TRUE,
                       verbose=TRUE)
  }

  #Apply to ASFRs
  asfr_sum <- asfr_sample_all <- data.frame()
  for (ctry in ctrycodes) {
    sample_pasfr <- read_tsv(paste0(tfr.dir,"/pasfr_traj/pasfr_",ctry,".txt")) %>%
      filter(age %in% agerange)
    tfrmat <- get.tfr.trajectories(tfr.pred,ctry)[-c(1),]
    sample_asfr <- sample_pasfr
    for (i in 1:500)
      sample_asfr[,paste(i)] <- rep(tfrmat[,i],each=30)*sample_pasfr[,paste(i)]/100
    sample_asfr_coh <- sample_asfr %>% mutate(coh=year-age)
    sample_asfr_sum <- sample_asfr %>%
      pivot_longer(-c("country_code","name","age","year"),
                   names_to="iter",values_to="ASFR")
    sample_asfr_sum %<>% group_by(name,age,year) %>%
      summarise(med=median(ASFR),
                lwr=quantile(ASFR,p=0.025),
                upr=quantile(ASFR,p=0.975)) %>%
      mutate(ctrycode=ctry)
    asfr_sum <- rbind(asfr_sum,sample_asfr_sum)
    asfr_sample_all <- rbind(asfr_sample_all,sample_asfr_coh)
  }
  write_tsv(asfr_sample_all,file=paste0("output/bayesTFR_asfr_sample_all_",joy,".txt"))
  
  #Plot ASFR prediction intervals for a given country
  for (i in 1:length(ctrycodes)) {
    png(file=paste0("plots/bayesTFR_ASFR_int_",ctrynames[i],"_",joy,".png"),width=30,height=15,units="cm",res=400,family="Calibri")
    print({asfr_sum %>% filter(ctrycode==ctrycodes[i]) %>%
        ggplot(aes(x=age,y=med,color=factor(year))) +
        coord_cartesian(ylim=c(0,0.2)) +
        labs(x="Age",y="ASFR") +
        scale_x_continuous(breaks=seq(15,44,5),minor_breaks=0:100) +
        geom_ribbon(aes(ymin=lwr,ymax=upr,fill=factor(year)),alpha=0.3, show.legend = F, linetype=2) +
        geom_line(aes(color=factor(year),group=factor(year)),show.legend = F) +
        facet_wrap(~year) + theme_bw() })
    dev.off()
  }
  
  #Plot ASFR posterior medians for all countries
  png(file=paste0("plots/bayesTFR_ASFR_med_",joy,".png"),width=30,height=15,units="cm",res=400,family="Calibri")
  print({asfr_sum %>% 
      ggplot(aes(x=age,y=med,color=year)) +
      scale_x_continuous(breaks=seq(15,44,5),minor_breaks=0:100) +
      scale_color_gradientn(colours=viridis(100)[100:1], guide = guide_colorbar(barheight = 10,frame.colour="black",ticks.colour="black"), breaks=seq(1900,2100,5),limits=c(joy+1,joy+30)) +
      geom_line(aes(color=year,group=year),show.legend = T,linewidth=0.3) +
      labs(x="Age",y="ASFR",color="Year",fill="Year") +
      coord_cartesian(ylim=c(0,0.2)) +
      facet_wrap(~name,nrow=5) +
      theme_bw() })
  dev.off()
}
