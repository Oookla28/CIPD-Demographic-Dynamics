###Assess exposure approximation method
##Read in rates
load(paste0("data/",ctry,"_rates.RData"))


##Extract register-based exposure estimates from HFD ('exposRRpac.txt' from 'parityexp.zip')
filenames <- as.character(unzip("data/parityexp.zip", list = TRUE)$Name)
unzip("data/parityexp.zip", filenames[which(filenames=="Files/zip_w/exposRRpac.txt")], exdir="data", junkpaths=T)


##Process cumulation-based exposure estimates
#For a given cohort c and age a, we assume that half the exposure is from
#year c+a (the year that cohort c turns age a), and half is from the
#year (c+a)+1 (the year that cohort c turns age a+1)

#Create copy of alldata for parities 0-3+ and selected columns
alldata2 <- alldata %>% filter(par>0) %>% select(par,Cohort,Age,N)
#Add a year column (c+a, the earlier year)
alldata2 %<>% mutate(Year=Cohort+Age)
#Add cohort quality issue indicator
if(!is.null(coh_iss[[ctry]])) alldata2 %<>% mutate(iss=ifelse(Cohort %in% coh_iss[[ctry]],T,F))
if(is.null(coh_iss[[ctry]])) alldata2$iss <- F
#Duplicate alldata2, with the duplication having the later year ((c+a)+1)
alldata2 <- rbind(alldata2,alldata2 %>% mutate(Year=Year+1))
#Halve the exposures
alldata2 %<>% mutate(N=0.5*N)
#Add a column to count the number of copies of each age-cohort combination
alldata2 %<>% group_by(par,Age,Year) %>% mutate(ind=n())
#Remove age-cohort combinations that are only partially observed
alldata2 %<>% filter(ind==2)
#Sum exposures and cohort quality issue indicator across parity, age and period
alldata2 %<>% group_by(par,Age,Year) %>% summarise(N=sum(N), iss=sum(as.numeric(iss)))
#Recode cohort quality issue variable
alldata2 %<>% mutate(iss=ifelse(iss==0,F,T))
#Arrange by year
alldata2 %<>% arrange(par,Year,Age)


##Process register-based exposure estimates
#Read in register-based exposure estimates
expos  <-  read.table("data/exposRRpac.txt",skip=2,as.is=T, header=T) %>% as_tibble()
#Extract countries in expos and our sample
inc_ctry <- intersect(ctrys,unique(expos$Code))

if (ctry %in% inc_ctry) {
  #Restrict to country and reproductive age range
  expos %<>% filter(Code == ctry, Age %in% 15:44)
  #Make into a long data frame
  expos %<>% pivot_longer(-c("Code","Year","Age"),names_to="par",values_to="N_reg")
  #Make par a numeric variable
  expos %<>% mutate(par=as.numeric(substr(par,2,2))+1)
  #Combine parities 3 and above
  expos %<>% mutate(par=ifelse(par>4,4,par))
  expos %<>% group_by(par,Year,Age) %>% summarise(N_reg=sum(N_reg))
  
  
  ##Merge datasets
  compdat <- alldata2 %>% left_join(expos)
  #Remove rows with missing exposures
  compdat %<>% filter(!is.na(N) & !is.na(N_reg))
  
  
  ##Extract age and year correlations
  compdat_age_cor_update <- compdat %>% group_by(par,Age) %>% summarise(cor=cor(N,N_reg)) %>% mutate(Code=ctry)
  compdat_yr_cor_update <- compdat %>% group_by(par,Year) %>% summarise(cor=cor(N,N_reg)) %>% mutate(Code=ctry)
  if (match(ctry,inc_ctry) == 1) {
    compdat_age_cor <- compdat_age_cor_update
    compdat_yr_cor <- compdat_yr_cor_update
  }
  if (match(ctry,inc_ctry) > 1) {
    compdat_age_cor <- rbind(compdat_age_cor,compdat_age_cor_update)
    compdat_yr_cor <- rbind(compdat_yr_cor,compdat_yr_cor_update)
  }
  
  ##Extract exposure estimates 
  if (match(ctry,inc_ctry) == 1) compdat_all <- compdat %>% mutate(Code=ctry)
  if (match(ctry,inc_ctry) > 1) compdat_all <- rbind(compdat_all,compdat %>% mutate(Code=ctry))
  
  
  ##Plot exposure estimates
  supp.labs1 <- c("ASFR","Parity 0","Parity 1","Parity 2","Parity 3+")
  names(supp.labs1) <- c(0,1,2,3,4)
  
  #By age
  png(file=paste0("plots/",ctry,"_comp_expos_age.png"),width=20,height=6.5,units="cm",res=400,family="Calibri")
  print({compdat %>% ggplot(aes(x=N/1000,y=N_reg/1000,color=Age)) +
      geom_point(data=compdat %>% filter(iss==T),size=0.75,shape=1) +
      geom_point(data=compdat %>% filter(iss==F),size=0.75,shape=19) +
      facet_wrap(~par,scales="free",nrow=1,labeller=labeller(par=supp.labs1)) +
      geom_abline(aes(intercept=0,slope=1)) +
      labs(x="Cumulation-based estimate (thousands)",y="Register-based estimate (thousands)") +
      scale_color_gradientn(colours=viridis(100)[100:1], guide = guide_colorbar(barheight = 10,frame.colour="black",ticks.colour="black"), breaks = seq(15,45,5), limits=c(15,44)) +
      theme_bw()})
  dev.off()
  
  #By year
  years <- unique(compdat$Year)
  nyears <- length(years)
  if (nyears<=7) break_values <- years
  if (nyears>7)  break_values <- seq(1800,2100,5)
  png(file=paste0("plots/",ctry,"_comp_expos_yr.png"),width=20,height=6.5,units="cm",res=400,family="Calibri")
  print({compdat %>% ggplot(aes(x=N/1000,y=N_reg/1000,color=Year)) +
      geom_point(data=compdat %>% filter(iss==T),size=0.75,shape=1) +
      geom_point(data=compdat %>% filter(iss==F),size=0.75,shape=19) +
      facet_wrap(~par,scales="free",nrow=1,labeller=labeller(par=supp.labs1)) +
      geom_abline(aes(intercept=0,slope=1)) +
      labs(x="Cumulation-based estimate (thousands)",y="Register-based estimate (thousands)") +
      scale_color_gradientn(colours=viridis(100)[100:1], guide = guide_colorbar(barheight = 10,frame.colour="black",ticks.colour="black"), breaks = break_values, limits=range(compdat$Year)) +
      theme_bw()})
  dev.off()
  
  ##Plots for all countries
  if (match(ctry,inc_ctry) == length(inc_ctry)) {
    #Plot exposure estimates for all countries
    #By age
    png(file=paste0("plots/ALL_comp_expos_age.png"),width=20,height=6.5,units="cm",res=400,family="Calibri")
    print({compdat_all %>% ggplot(aes(x=N/1000,y=N_reg/1000,color=Age)) +
        geom_point(data=compdat_all %>% filter(iss==T),size=0.75,shape=1) +
        geom_point(data=compdat_all %>% filter(iss==F),size=0.75,shape=19) +
        facet_wrap(~par,scales="free",nrow=1,labeller=labeller(par=supp.labs1)) +
        geom_abline(aes(intercept=0,slope=1)) +
        labs(x="Cumulation-based estimate (thousands)",y="Register-based estimate (thousands)") +
        scale_color_gradientn(colours=viridis(100)[100:1], guide = guide_colorbar(barheight = 10,frame.colour="black",ticks.colour="black"), breaks = seq(15,45,5), limits=c(15,44)) +
        theme_bw()})
    dev.off()
    
    #By year
    years <- unique(compdat_all$Year)
    nyears <- length(years)
    if (nyears<=7) break_values <- years
    if (nyears>7)  break_values <- seq(1800,2100,5)
    png(file=paste0("plots/ALL_comp_expos_yr.png"),width=20,height=6.5,units="cm",res=400,family="Calibri")
    print({compdat_all %>% ggplot(aes(x=N/1000,y=N_reg/1000,color=Year)) +
        geom_point(data=compdat_all %>% filter(iss==T),size=0.75,shape=1) +
        geom_point(data=compdat_all %>% filter(iss==F),size=0.75,shape=19) +
        facet_wrap(~par,scales="free",nrow=1,labeller=labeller(par=supp.labs1)) +
        geom_abline(aes(intercept=0,slope=1)) +
        labs(x="Cumulation-based estimate (thousands)",y="Register-based estimate (thousands)") +
        scale_color_gradientn(colours=viridis(100)[100:1], guide = guide_colorbar(barheight = 10,frame.colour="black",ticks.colour="black"), breaks = break_values, limits=range(compdat_all$Year)) +
        theme_bw()})
    dev.off()
    
    #By country
    png(file=paste0("plots/ALL_comp_expos_ctry.png"),width=20,height=10,units="cm",res=400,family="Calibri")
    print({compdat_all %>% mutate(ctry=ctrys_long[Code]) %>%
        ggplot(aes(x=N/1000,y=N_reg/1000,color=ctry)) +
        geom_point(data=compdat_all %>% mutate(ctry=ctrys_long[Code]) %>% filter(iss==T),size=0.75,shape=1) +
        geom_point(data=compdat_all %>% mutate(ctry=ctrys_long[Code]) %>% filter(iss==F),size=0.75,shape=19) +
        facet_wrap(~par,scales="free",nrow=1,labeller=labeller(par=supp.labs1)) +
        geom_abline(aes(intercept=0,slope=1)) +
        guides(color=guide_legend(nrow=2,byrow=T)) +
        labs(x="Cumulation-based estimate (thousands)",y="Register-based estimate (thousands)") +
        theme_bw() +
        theme(legend.position="bottom")})
    dev.off()
    
    #Plot correlations for all countries
    #By age
    png(file=paste0("plots/ALL_comp_expos_cor_age.png"),width=20,height=10,units="cm",res=400,family="Calibri")
    print({compdat_age_cor %>% mutate(ctry=ctrys_long[Code]) %>%
        ggplot(aes(x=Age,y=cor,color=ctry)) +
        geom_point(size=0.75) +
        facet_grid(~par,labeller=labeller(par=supp.labs1)) +
        scale_x_continuous(breaks=seq(15,40,5),minor_breaks=0:100) +
        labs(x="Age",y="Correlation",color="Country") +
        guides(color=guide_legend(nrow=2,byrow=T)) +
        theme_bw() +
        theme(legend.position="bottom")})
    dev.off()
    
    #By year
    png(file=paste0("plots/ALL_comp_expos_cor_yr.png"),width=20,height=10,units="cm",res=400,family="Calibri")
    print({compdat_yr_cor %>% mutate(ctry=ctrys_long[Code]) %>%
        ggplot(aes(x=Year,y=cor,color=ctry)) +
        geom_point(size=0.75) +
        facet_wrap(~par,scales="free",nrow=1,labeller=labeller(par=supp.labs1)) +
        scale_x_continuous(breaks=seq(1900,2100,20),minor_breaks=1900:2100) +
        labs(x="Year",y="Correlation",color="Country") +
        guides(color=guide_legend(nrow=2,byrow=T)) +
        theme_bw() +
        theme(legend.position="bottom")})
    dev.off()
  }
}
