###Plot model output for multiple countries/JOYs
##Set-up
supp.lab.p <- c("ASFR",paste0("Parity ",c(0,1,2,"3+"))); names(supp.lab.p) <- 0:4
supp.lab.ppr <- c("PPR 0 -> 1","PPR 1 -> 2","PPR 2 -> 3","PPR 3+ -> 4+"); names(supp.lab.ppr) <- 2:5


##Plot fitted rates with observed rates for two countries and one JOY
ctry2 <- c("JPN","LTU") #Select two countries
#ctry2 <- c("DNK","IRL") #Select two countries
joy <- joys[4] #Select a JOY
ctryvec <- rep(ctry2,2)
ctrylvec <- rep(ctrys_long[ctry2],2)
mods <- c("Aggregate","Parity-specific")
modvec <- rep(mods,each=length(ctry2))
files <- paste0("mod_",rep(c("a","p"),each=length(ctry2)),"_",ctryvec,"_",joy,".RData")
redvec <- c(red_sel_lst[[ctry2[1]]][which(joy_sel_lst[[ctry2[1]]]==joy)],
            red_sel_lst[[ctry2[2]]][which(joy_sel_lst[[ctry2[2]]]==joy)])

for (ctry in ctry2) {
  red <- redvec[which(ctry2==ctry)]
  load(paste0("data/",ctry,"_rates.RData"))
  load(paste0("output/",files[which(ctryvec==ctry)[1]])); stanouta <- stanout
  load(paste0("output/",files[which(ctryvec==ctry)[2]])); stanoutp <- stanout
  rate <- cbind(extract(stanouta,"asfr",permuted=F)[,1,],
                extract(stanoutp,pars=c("asfr",paste0("pf",1:4)),permuted=F)[,1,])
  if (red==F) cohrange <- (joy-max(agerange)-10):(joy-min(agerange))
  if (red==T) cohrange <- (joy-max(agerange)-5):(joy-min(agerange))
  Nf <- length(agerange)*length(cohrange)
  ratedat <- data.frame(ctry_l=ctrys_long[ctry][[1]],
                        expand.grid(Age=agerange,Cohort=cohrange),
                        mod=rep(mods,c(Nf,5*Nf)),
                        par=rep(c(0,0,1:4),each=Nf),
                        med=apply(rate,2,median))
  ratedat %<>% left_join(alldata %>% mutate(par=as.numeric(par)) %>% select(par,Cohort,Age,x))
  ratedat <- rbind(ratedat,ratedat) %>%
    mutate(type=rep(c("Fitted + forecast","Observed"),each=n()/2))
  ratedat$med[ratedat$type=="Observed"] <- ratedat$x[ratedat$type=="Observed"]
  ratedatf <- ratedatp <- ratedat
  ratedatp$med[which(ratedatp$Age+ratedatp$Cohort>joy)] <- NA
  ratedatf$med[which(ratedatf$Age+ratedatf$Cohort<joy)] <- NA
  if (ctry==ctry2[1]) ratedat_ctry <- ratedat else ratedat_ctry <- rbind(ratedat_ctry,ratedat)
  if (ctry==ctry2[1]) ratedatp_ctry <- ratedatp else ratedatp_ctry <- rbind(ratedatp_ctry,ratedatp)
  if (ctry==ctry2[1]) ratedatf_ctry <- ratedatf else ratedatf_ctry <- rbind(ratedatf_ctry,ratedatf)
}

png(file=paste0("plots/mod_ap_",ctry2[1],"_",ctry2[2],"_rates_med_",joy,"_2.png"),width=20,height=15.75,units="cm",res=400,family="Calibri")
print({ ggplot(ratedatp_ctry %>% mutate(ctry_l=factor(ctry_l,levels=ctrys_long[ctry2])) %>%
                 filter(med<1),aes(x=Age,color=Cohort)) +
    facet_nested(ctry_l+type~mod+par,labeller=labeller(par=supp.lab.p)) +
    geom_line(aes(y=med,group=Cohort),linewidth=0.2) +
    geom_line(data=ratedatf_ctry %>% mutate(ctry_l=factor(ctry_l,levels=ctrys_long[ctry2])) %>%
                filter(med<1),aes(y=med,group=Cohort),linewidth=0.2,linetype=5) +
    coord_cartesian(ylim=c(0,0.3)) +
    scale_x_continuous(breaks=seq(15,45,5),minor_breaks=agerange,expand=c(0.025,0.025)) +
    labs(x="Age",y="Rate",color="Cohort") +
    scale_color_gradientn(colours=viridis(100)[100:1], guide = guide_colorbar(barwidth=20, barheight = 1,frame.colour="black",ticks.colour="black"), breaks = seq(1920,2000,5), limits=range(cohrange)) +
    theme_bw() + theme(text = element_text("Calibri"),axis.text=element_text(color="black"),
                       axis.title=element_text(face="bold"),legend.position = "bottom") })
dev.off()

tiff(file=paste0("plots/mod_ap_",ctry2[1],"_",ctry2[2],"_rates_med_",joy,"_2.tif"),
     width=20,height=15.75,units="cm",res=600,family="Calibri",compression="none")
print({ ggplot(ratedatp_ctry %>% mutate(ctry_l=factor(ctry_l,levels=ctrys_long[ctry2])) %>%
                 filter(med<1),aes(x=Age,color=Cohort)) +
    facet_nested(ctry_l+type~mod+par,labeller=labeller(par=supp.lab.p)) +
    geom_line(aes(y=med,group=Cohort),linewidth=0.2) +
    geom_line(data=ratedatf_ctry %>% mutate(ctry_l=factor(ctry_l,levels=ctrys_long[ctry2])) %>%
                filter(med<1),aes(y=med,group=Cohort),linewidth=0.2,linetype=5) +
    coord_cartesian(ylim=c(0,0.3)) +
    scale_x_continuous(breaks=seq(15,45,5),minor_breaks=agerange,expand=c(0.025,0.025)) +
    labs(x="Age",y="Rate",color="Cohort") +
    scale_color_gradientn(colours=viridis(100)[100:1], guide = guide_colorbar(barwidth=20, barheight = 1,frame.colour="black",ticks.colour="black"), breaks = seq(1920,2000,5), limits=range(cohrange)) +
    theme_bw() + theme(text = element_text("Calibri"),axis.text=element_text(color="black"),
                       axis.title=element_text(face="bold"),legend.position = "bottom") })
dev.off()


##Generate uncertainty for four countries and one JOY
ctry2 <- c("DNK","JPN","IRL","LTU") #Select four countries
joy <- joys[4] #Select a JOY
ctryvec <- rep(ctry2,2)
ctrylvec <- rep(ctrys_long[ctry2],2)
mods <- c("Aggregate","Parity-specific")
modvec <- rep(mods,each=length(ctry2))
files <- paste0("mod_",rep(c("a","p"),each=length(ctry2)),"_",ctryvec,"_",joy,".RData")
redvec <- rep(c(red_sel_lst[[ctry2[1]]][which(joy_sel_lst[[ctry2[1]]]==joy)],
            red_sel_lst[[ctry2[2]]][which(joy_sel_lst[[ctry2[2]]]==joy)],
            red_sel_lst[[ctry2[3]]][which(joy_sel_lst[[ctry2[3]]]==joy)],
            red_sel_lst[[ctry2[4]]][which(joy_sel_lst[[ctry2[4]]]==joy)]),2)

genunc <- gen_unc(files,agerange,rep(joy,length(files)),modvec,ctryvec,redvec)
moddat <- data.frame(mod=factor(1:length(files)),ctry=ctryvec,ctryl=ctrylvec,type=modvec)
intdat <- data.frame(joy2=seq(2000,2015,5),int=seq(2000,2015,5)-max(agerange))


##Combine datasets
alldata2 <- data.frame()
for (i in ctry2) {
  load(paste0("data/",i,"_rates.RData"))
  alldata2 <- rbind(alldata2,data.frame(ctry=i,alldata))
}

#Compute observed values of summary measures
cfrobs <- alldata2 %>% filter(par==0) %>%
  group_by(ctry,Cohort) %>% mutate(mina=min(Age),maxa=max(Age)) %>%
  filter(mina==min(agerange),maxa==max(agerange)) %>%
  summarise(obs=sum(x)) %>% rename("coh"=Cohort)
pprobs <- alldata2 %>% filter(par>0) %>%
  left_join(alldata2 %>% filter(par==0) %>% select(ctry,Cohort,Age,N) %>% rename("Ntot"=N), multiple="all") %>%
  mutate(uasfr=n/Ntot) %>% select(-N,-n,-x,-Ntot) %>%
  pivot_wider(names_from=par,values_from="uasfr",names_prefix="asfr") %>%
  group_by(ctry,Cohort) %>%
  summarise(mina=min(Age),maxa=max(Age),
            cfr1=sum(asfr1,na.rm=T),cfr2=sum(asfr2,na.rm=T),
            cfr3=sum(asfr3,na.rm=T),cfr4=sum(asfr4,na.rm=T),
            ppr1=cfr1,ppr2=cfr2/cfr1,ppr3=cfr3/cfr2,ppr4=cfr4/(cfr3+cfr4)) %>%
  select(ctry,Cohort,starts_with("ppr"),mina,maxa) %>%
  pivot_longer(-c("ctry","Cohort","mina","maxa")) %>%
  filter(mina==15,maxa==44) %>% select(-mina,-maxa) %>%
  rename("ind"=name,"obs"=value,"coh"=Cohort) %>%
  mutate(ind=as.numeric(substr(ind,4,4))+1)
fsdobs <- alldata2 %>% filter(par>0) %>%
  left_join(alldata2 %>% filter(par==0) %>% select(ctry,Cohort,Age,N) %>% rename("Ntot"=N), multiple="all") %>%
  mutate(uasfr=n/Ntot) %>% select(-N,-n,-x,-Ntot) %>%
  pivot_wider(names_from=par,values_from="uasfr",names_prefix="asfr") %>%
  group_by(ctry,Cohort) %>%
  summarise(mina=min(Age),maxa=max(Age),
            cfr1=sum(asfr1,na.rm=T),cfr2=sum(asfr2,na.rm=T),
            cfr3=sum(asfr3,na.rm=T),cfr4=sum(asfr4,na.rm=T),
            ppr1=cfr1,ppr2=cfr2/cfr1,ppr3=cfr3/cfr2,ppr4=cfr4/(cfr3+cfr4),
            fsd1=1-cfr1,fsd2=cfr1-cfr2,fsd3=cfr2-cfr3,fsd4=cfr3) %>%
  select(ctry,Cohort,starts_with("fsd"),mina,maxa) %>%
  pivot_longer(-c("ctry","Cohort","mina","maxa")) %>%
  filter(mina==15,maxa==44) %>% select(-mina,-maxa) %>%
  rename("ind"=name,"obs"=value,"coh"=Cohort) %>%
  mutate(ind=as.numeric(substr(ind,4,4))+1)


##Plot CFR posterior distribution
png(file=paste0("plots/mod_ap_cfr_post_",joy,".png"),width=15,height=6.75,units="cm",res=400,family="Calibri")
print({genunc$cfr %>% left_join(moddat) %>% left_join(cfrobs,multiple="all") %>%
    ggplot(aes(x=coh,color=type,fill=type)) +
    geom_hline(aes(yintercept=2),color="grey50") +
    geom_ribbon(aes(ymin=lwr,ymax=upr),linetype=2,linewidth=0.3,alpha=0.3) +
    geom_line(aes(y=med)) +
    facet_grid(~ctryl,scales="free_x") +
    scale_x_continuous(breaks=seq(1900,2005,5),minor_breaks=1900:2100,expand=c(0.02,0.02)) +
    scale_y_continuous(breaks=seq(0,3,0.2),limits=c(1.25,2.4)) +
    geom_point(aes(x=coh,y=obs,fill=type),color=1,size=0.5,show.legend=F) +
    geom_vline(data=intdat %>% filter(joy2==joy),aes(xintercept=int)) +
    labs(x="Cohort",y="CFR",fill="Model",color="Model") +
    theme_bw() + theme(text = element_text("Calibri"),legend.position = "bottom",
                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) })
dev.off()


##Plot PPR posterior distribution
png(file=paste0("plots/mod_ap_ppr_post_",joy,".png"),width=15,height=10,units="cm",res=400,family="Calibri")
print({genunc$ppr %>% left_join(moddat) %>% left_join(pprobs,multiple="all") %>%
    left_join(data.frame(ind=2:5,ppr=c("PPR 0 -> 1","PPR 1 -> 2","PPR 2 -> 3","PPR 3+ -> 4+"))) %>%
    ggplot(aes(x=coh,color=factor(ppr),fill=factor(ppr))) +
    geom_ribbon(aes(ymin=lwr,ymax=upr,group=ppr),linetype=2,linewidth=0.3,alpha=0.3) +
    geom_line(aes(y=med,group=ppr)) +
    facet_grid(~ctryl,scales="free") +
    coord_cartesian(ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(1900,2005,5),minor_breaks=1900:2100,expand=c(0.02,0.02)) +
    scale_y_continuous(breaks=seq(0,3,0.1)) +
    geom_point(aes(x=coh,y=obs,color=factor(ppr)),pch=21,fill="black",stroke=0.5,size=0.75,show.legend=F) +
    geom_vline(data=intdat %>% filter(joy2==joy),aes(xintercept=int)) +
    labs(x="Cohort",y="Posterior distribution",fill="Measure",color="Measure") +
    theme_bw() + theme(text = element_text("Calibri"),legend.position = "bottom",
                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))})
dev.off()


##Plot CFR and PPR posterior distributions together
intdat2 <- data.frame(type2="CFR",yintercept=2)
png(file=paste0("plots/mod_ap_cfr_ppr_post_",joy,".png"),width=15,height=14,units="cm",res=400,family="Calibri")
print({rbind(genunc$cfr %>% left_join(moddat) %>% left_join(cfrobs,multiple="all") %>%
               mutate(type2="CFR",ind=ifelse(type=="Aggregate",0,1)),
             genunc$ppr %>% left_join(moddat) %>% left_join(pprobs,multiple="all") %>% mutate(type2="PPRs")) %>%
    left_join(data.frame(ind=0:5,ppr=c("CFR - Aggregate model","CFR - Parity-specific model","PPR 0 -> 1","PPR 1 -> 2","PPR 2 -> 3","PPR 3+ -> 4+"))) %>%
    ggplot(aes(x=coh,color=factor(ppr),fill=factor(ppr))) +
    geom_hline(data=intdat2,aes(yintercept=yintercept),color="grey50") +
    geom_ribbon(aes(ymin=lwr,ymax=upr,group=ppr),linetype=2,linewidth=0.3,alpha=0.3) +
    geom_line(aes(y=med,group=ppr)) +
    facet_grid(type2~ctryl,scales="free") +
    force_panelsizes(rows = c(1,2)) +
    scale_fill_manual(values=hue_pal()(6)[c(1,4,2,5,6,3)]) +
    scale_color_manual(values=hue_pal()(6)[c(1,4,2,5,6,3)]) +
    scale_x_continuous(breaks=seq(1900,2005,5),minor_breaks=1900:2100,expand=c(0.02,0.02)) +
    scale_y_continuous(breaks=c(seq(0,1,0.1),seq(1.2,3,0.2))) +
    geom_point(aes(x=coh,y=obs,color=factor(ppr)),pch=21,fill="black",stroke=0.5,size=0.75,show.legend=F) +
    geom_vline(data=intdat %>% filter(joy2==joy),aes(xintercept=int)) +
    labs(x="Cohort",y="Posterior Distribution",fill="Measure",color="Measure") +
    theme_bw() + theme(text = element_text("Calibri"),legend.position = "bottom",
                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                       axis.text=element_text(color="black"),
                       axis.title=element_text(face="bold"))})
dev.off()

tiff(file=paste0("plots/mod_ap_cfr_ppr_post_",joy,".tif"),width=15,height=14,units="cm",res=600,family="Calibri",compression="none")
print({rbind(genunc$cfr %>% left_join(moddat) %>% left_join(cfrobs,multiple="all") %>%
               mutate(type2="CFR",ind=ifelse(type=="Aggregate",0,1)),
             genunc$ppr %>% left_join(moddat) %>% left_join(pprobs,multiple="all") %>% mutate(type2="PPRs")) %>%
    left_join(data.frame(ind=0:5,ppr=c("CFR - Aggregate model","CFR - Parity-specific model","PPR 0 -> 1","PPR 1 -> 2","PPR 2 -> 3","PPR 3+ -> 4+"))) %>%
    ggplot(aes(x=coh,color=factor(ppr),fill=factor(ppr))) +
    geom_hline(data=intdat2,aes(yintercept=yintercept),color="grey50") +
    geom_ribbon(aes(ymin=lwr,ymax=upr,group=ppr),linetype=2,linewidth=0.3,alpha=0.3) +
    geom_line(aes(y=med,group=ppr)) +
    facet_grid(type2~ctryl,scales="free") +
    force_panelsizes(rows = c(1,2)) +
    scale_fill_manual(values=hue_pal()(6)[c(1,4,2,5,6,3)]) +
    scale_color_manual(values=hue_pal()(6)[c(1,4,2,5,6,3)]) +
    scale_x_continuous(breaks=seq(1900,2005,5),minor_breaks=1900:2100,expand=c(0.02,0.02)) +
    scale_y_continuous(breaks=c(seq(0,1,0.1),seq(1.2,3,0.2))) +
    geom_point(aes(x=coh,y=obs,color=factor(ppr)),pch=21,fill="black",stroke=0.5,size=0.75,show.legend=F) +
    geom_vline(data=intdat %>% filter(joy2==joy),aes(xintercept=int)) +
    labs(x="Cohort",y="Posterior Distribution",fill="Measure",color="Measure") +
    theme_bw() + theme(text = element_text("Calibri"),legend.position = "bottom",
                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                       axis.text=element_text(color="black"),
                       axis.title=element_text(face="bold"))})
dev.off()


##Plot CFR, PPR and FSD posterior distributions together
intdat2 <- data.frame(type2="CFR",yintercept=2); intdat2$type2 <- factor(intdat2$type2,levels=c("CFR","PPRs","FSD"))
png(file=paste0("plots/mod_ap_cfr_ppr_fsd_post_",joy,".png"),width=15,height=19,units="cm",res=400,family="Calibri")
print({rbind(genunc$cfr %>% left_join(moddat) %>% left_join(cfrobs,multiple="all") %>%
               mutate(type2="CFR",ind=ifelse(type=="Aggregate",0,1)),
             genunc$ppr %>% left_join(moddat) %>% left_join(pprobs,multiple="all") %>% mutate(type2="PPRs"),
             genunc$fsd %>% left_join(moddat) %>% left_join(fsdobs,multiple="all") %>% mutate(type2="FSD")) %>%
    left_join(data.frame(ind=0:5,ppr=c("CFR - Aggregate model","CFR - Parity-specific model",
                                       "PPR 0 -> 1\nFSD: 0 children",
                                       "PPR 1 -> 2\nFSD: 1 child",
                                       "PPR 2 -> 3\nFSD: 2 children",
                                       "PPR 3+ -> 4+\nFSD: 3+ children"))) %>%
    mutate(type2=factor(type2,levels=c("CFR","PPRs","FSD"))) %>%
    ggplot(aes(x=coh,color=factor(ppr),fill=factor(ppr))) +
    geom_hline(data=intdat2,aes(yintercept=yintercept),color="grey50") +
    geom_ribbon(aes(ymin=lwr,ymax=upr,group=ppr),linetype=2,linewidth=0.3,alpha=0.3) +
    geom_line(aes(y=med,group=ppr)) +
    facet_grid(type2~ctryl,scales="free") +
    force_panelsizes(rows = c(1,2,2)) +
    scale_fill_manual(values=hue_pal()(6)[c(1,4,2,5,6,3)]) +
    scale_color_manual(values=hue_pal()(6)[c(1,4,2,5,6,3)]) +
    scale_x_continuous(breaks=seq(1900,2005,5),minor_breaks=1900:2100,expand=c(0.02,0.02)) +
    scale_y_continuous(breaks=c(seq(0,1,0.1),seq(1.2,3,0.2))) +
    geom_point(aes(x=coh,y=obs,color=factor(ppr)),pch=21,fill="black",stroke=0.5,size=0.75,show.legend=F) +
    geom_vline(data=intdat %>% filter(joy2==joy),aes(xintercept=int)) +
    labs(x="Cohort",y="Posterior distribution",fill="Measure",color="Measure") +
    theme_bw() + theme(text = element_text("Calibri"),legend.position = "bottom",
                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))})
dev.off()
