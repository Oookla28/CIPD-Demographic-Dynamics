###Assess exposure approximation method
##Read in rates and combine
alldata2 <- data.frame()
for (ctry in ctrys) {
  load(paste0("data/",ctry,"_rates.RData"))
  alldata2 <- rbind(alldata2,alldata %>% mutate(ctry=ctry))
}

#Restrict to aggregate data and cohort range of interest
alldata2 %<>% filter(par==0,Cohort %in% 1946:2000)


##Compute ratio of exposures for a given cohort and country, and consecutive ages
alldata2 %<>% group_by(ctry,Cohort) %>%
  mutate(Nlag=lag(N),
         Nrat=N/Nlag) %>%
  filter(!is.na(Nrat))


##Plot exposure ratios
#Density plot
png(file=paste0("plots/ALL_comp_expos_scaling_dens.png"),width=15,height=7.5,units="cm",res=400,family="Calibri")
print({alldata2 %>% 
  ggplot(aes(x=Nrat)) +
  geom_density() +
  labs(x="Ratio of exposures",y="Density") +
  geom_vline(aes(xintercept=1),linetype=2) +
  theme_bw()})
dev.off()

#By age
png(file=paste0("plots/ALL_comp_expos_scaling_age.png"),width=20,height=10,units="cm",res=400,family="Calibri")
print({alldata2 %>% ggplot(aes(x=N/1000,y=Nlag/1000,color=Age)) +
    geom_point(size=0.75,shape=19) +
    geom_abline(aes(intercept=0,slope=1)) +
    labs(x="Current age exposure (thousands)",y="Previous age exposure (thousands)") +
    scale_color_gradientn(colours=viridis(100)[100:1], guide = guide_colorbar(barheight = 10,frame.colour="black",ticks.colour="black"), breaks = seq(15,45,5), limits=c(15,44)) +
    theme_bw()})
dev.off()

#By cohort
png(file=paste0("plots/ALL_comp_expos_scaling_yr.png"),width=20,height=10,units="cm",res=400,family="Calibri")
print({alldata2 %>% ggplot(aes(x=N/1000,y=Nlag/1000,color=Cohort)) +
    geom_point(size=0.75,shape=19) +
    geom_abline(aes(intercept=0,slope=1)) +
    labs(x="Current age exposure (thousands)",y="Previous age exposure (thousands)") +
    scale_color_gradientn(colours=viridis(100)[100:1], guide = guide_colorbar(barheight = 10,frame.colour="black",ticks.colour="black"), limits=range(alldata2$Cohort)) +
    theme_bw()})
dev.off()

#By country
png(file=paste0("plots/ALL_comp_expos_scaling_ctry.png"),width=30,height=15,units="cm",res=400,family="Calibri")
print({alldata2 %>% ggplot(aes(x=N/1000,y=Nlag/1000,color=ctry)) +
    geom_point(size=0.75,shape=19,show.legend = F) +
    geom_abline(aes(intercept=0,slope=1),linewidth=0.3) +
    facet_wrap(~ctry,nrow=3) +
    labs(x="Current age exposure (thousands)",y="Previous age exposure (thousands)") +
    theme_bw()})
dev.off()
