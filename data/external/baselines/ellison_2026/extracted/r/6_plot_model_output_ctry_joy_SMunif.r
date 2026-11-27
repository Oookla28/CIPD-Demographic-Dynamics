###Plot model output for a single country and JOY
##Setup
ctryvec <- rep(ctry,2)
ctrylvec <- rep(ctrys_long[ctry],2)
mods <- modvec <- c("Aggregate","Parity-specific")
if (red==F) files <- paste0("mod_",c("a","p"),"_SMunif_",ctry,"_",joy,".RData")
if (red==T) files <- paste0("mod_",c("a","p"),"_SMunif_",ctry,"_",joy,"_red.RData")
supp.lab.p <- c("ASFR",paste0("Parity ",c(0,1,2,"3+"))); names(supp.lab.p) <- 0:4
supp.lab.ppr <- c("PPR 0 -> 1","PPR 1 -> 2","PPR 2 -> 3","PPR 3+ -> 4+"); names(supp.lab.ppr) <- 2:5


##Load stan data
if (red==F) load(paste0("data/standata_",ctry,"_",joy,".RData"))
if (red==T) load(paste0("data/standata_",ctry,"_",joy,"_red.RData"))


##Plot basis function coefficient estimates
load(paste0("output/",files[1])); stanouta <- stanout
load(paste0("output/",files[2])); stanoutp <- stanout
beta <- cbind(extract(stanouta,"beta0",permuted=F)[,1,],
              extract(stanout,pars=paste0("beta",1:4),permuted=F)[,1,])
NbC <- standata$NdC+1
NbA <- floor(standata$Na/5)+2
betadat <- data.frame(mod=rep(mods,c(standata$betadim,4*standata$betadim)),
                      type=rep(0:4,each=standata$betadim),
                      x=rep(rep(1:NbA,NbC),5),
                      y=rep(rep(1:NbC,each=NbA),5),
                      med=apply(beta,2,median))
plotname <- ifelse(red==F,paste0("plots/mod_ap_SMunif_",ctry,"_beta_med_",joy,".png"),
                   paste0("plots/mod_ap_SMunif_",ctry,"_beta_med_",joy,"_red.png"))
png(file=plotname,width=15,height=10,units="cm",res=400,family="Calibri")
print({ggplot(betadat,aes(x=y,y=x,fill=med)) +
    facet_wrap(~mod+type, nrow=2, labeller = labeller(type=supp.lab.p)) +
    geom_tile() +
    scale_x_continuous(expand=c(0,0),breaks=1:20) +
    scale_y_continuous(expand=c(0,0),breaks=1:8) +
    labs(x="Increasing cohort ->",y="Increasing age ->",fill="Value") +
    scale_fill_gradientn(colours=rainbow(100, start=0, end=1/6), guide = guide_colorbar(barheight = 15,frame.colour="black",ticks.colour="black"),breaks=seq(100,-25,-2.5)) +
    theme_bw() + theme(text = element_text("Calibri"))})
dev.off()


##Plot fitted rates with observed rates
load(paste0("data/",ctry,"_rates.RData"))
load(paste0("output/",files[1])); stanouta <- stanout
load(paste0("output/",files[2])); stanoutp <- stanout
rate <- cbind(extract(stanouta,"asfr",permuted=F)[,1,],
              extract(stanoutp,pars=c("asfr",paste0("pf",1:4)),permuted=F)[,1,])
Nf <- length(agerange)*length(cohrange)
ratedat <- data.frame(expand.grid(Age=agerange,Cohort=cohrange),
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

plotname <- ifelse(red==F,paste0("plots/mod_ap_SMunif_",ctry,"_rates_med_",joy,"_2.png"),
                   paste0("plots/mod_ap_SMunif_",ctry,"_rates_med_",joy,"_red_2.png"))
png(file=plotname,width=20,height=10,units="cm",res=400,family="Calibri")
print({ ggplot(ratedatp %>% filter(med<1),aes(x=Age,color=Cohort)) +
    facet_grid(type~mod+par,labeller=labeller(par=supp.lab.p)) +
    geom_line(aes(y=med,group=Cohort),linewidth=0.2) + # 0.3 previously
    geom_line(data=ratedatf %>% filter(med<1),aes(y=med,group=Cohort),linewidth=0.2,linetype=5) +
    coord_cartesian(ylim=c(0,0.3)) +
    scale_x_continuous(breaks=seq(15,45,5),minor_breaks=agerange,expand=c(0.025,0.025)) +
    labs(x="Age",y="Rate",color="Cohort") +
    scale_color_gradientn(colours=viridis(100)[100:1], guide = guide_colorbar(barwidth=20, barheight = 1,frame.colour="black",ticks.colour="black"), breaks = seq(1920,2100,5), limits=range(cohrange)) +
    theme_bw() + theme(text = element_text("Calibri"),legend.position = "bottom") })
dev.off()


##Generate uncertainty
genunc <- gen_unc(files,agerange,rep(joy,2),modvec,ctryvec)
moddat <- data.frame(mod=factor(1:length(files)),ctry=ctryvec,ctryl=ctrylvec,type=modvec)
intdat <- data.frame(joy2=joys,int=joys-max(agerange))


##Combine datasets
load(paste0("data/",ctry,"_rates.RData"))
alldata2 <- alldata

#Compute observed values of summary measures
cfrobs <- alldata2 %>% filter(par==0) %>%
  group_by(Cohort) %>% mutate(mina=min(Age),maxa=max(Age)) %>%
  filter(mina==min(agerange),maxa==max(agerange)) %>%
  summarise(obs=sum(x)) %>% rename("coh"=Cohort)
pprobs <- alldata2 %>% filter(par>0) %>%
  left_join(alldata2 %>% filter(par==0) %>% select(Cohort,Age,N) %>% rename("Ntot"=N), multiple="all") %>%
  mutate(uasfr=n/Ntot) %>% select(-N,-n,-x,-Ntot) %>%
  pivot_wider(names_from=par,values_from="uasfr",names_prefix="asfr") %>%
  group_by(Cohort) %>%
  summarise(mina=min(Age),maxa=max(Age),
            cfr1=sum(asfr1,na.rm=T),cfr2=sum(asfr2,na.rm=T),
            cfr3=sum(asfr3,na.rm=T),cfr4=sum(asfr4,na.rm=T),
            ppr1=cfr1,ppr2=cfr2/cfr1,ppr3=cfr3/cfr2,ppr4=cfr4/(cfr3+cfr4)) %>%
  select(Cohort,starts_with("ppr"),mina,maxa) %>%
  pivot_longer(-c("Cohort","mina","maxa")) %>%
  filter(mina==min(agerange),maxa==max(agerange)) %>% select(-mina,-maxa) %>%
  rename("ind"=name,"obs"=value,"coh"=Cohort) %>%
  mutate(ind=as.numeric(substr(ind,4,4))+1)
fsdobs <- alldata2 %>% filter(par>0) %>%
  left_join(alldata2 %>% filter(par==0) %>% select(Cohort,Age,N) %>% rename("Ntot"=N), multiple="all") %>%
  mutate(uasfr=n/Ntot) %>% select(-N,-n,-x,-Ntot) %>%
  pivot_wider(names_from=par,values_from="uasfr",names_prefix="asfr") %>%
  group_by(Cohort) %>%
  summarise(mina=min(Age),maxa=max(Age),
            cfr1=sum(asfr1,na.rm=T),cfr2=sum(asfr2,na.rm=T),
            cfr3=sum(asfr3,na.rm=T),cfr4=sum(asfr4,na.rm=T),
            ppr1=cfr1,ppr2=cfr2/cfr1,ppr3=cfr3/cfr2,ppr4=cfr4/(cfr3+cfr4),
            fsd1=1-cfr1,fsd2=cfr1-cfr2,fsd3=cfr2-cfr3,fsd4=cfr3) %>%
  select(Cohort,starts_with("fsd"),mina,maxa) %>%
  pivot_longer(-c("Cohort","mina","maxa")) %>%
  filter(mina==min(agerange),maxa==max(agerange)) %>% select(-mina,-maxa) %>%
  rename("ind"=name,"obs"=value,"coh"=Cohort) %>%
  mutate(ind=as.numeric(substr(ind,4,4))+1)


##Plot CFR posterior distribution
plotname <- ifelse(red==F,paste0("plots/mod_ap_SMunif_",ctry,"_cfr_post_",joy,".png"),
                   paste0("plots/mod_ap_SMunif_",ctry,"_cfr_post_",joy,"_red.png"))
png(file=plotname,width=10,height=10,units="cm",res=400,family="Calibri")
print({genunc$cfr %>% left_join(moddat) %>% left_join(cfrobs,multiple="all") %>%
    ggplot(aes(x=coh,color=type,fill=type)) +
    geom_hline(aes(yintercept=2),color="grey50") +
    geom_ribbon(aes(ymin=lwr,ymax=upr),linetype=2,linewidth=0.3,alpha=0.3) +
    geom_line(aes(y=med)) +
    facet_grid(~paste0(ctryl,", ",joy),scales="free_x") +
    scale_x_continuous(breaks=seq(1900,2005,5),minor_breaks=1900:2100,expand=c(0.02,0.02)) +
    coord_cartesian(ylim=c(1.25,2.3)) +
    scale_y_continuous(breaks=seq(0,3,0.2)) +
    geom_point(aes(x=coh,y=obs,fill=NULL),color=1,size=0.75,show.legend=F) +
    geom_vline(data=intdat %>% filter(joy2==joy),aes(xintercept=int)) +
    labs(x="Cohort",y="CFR",fill="Model",color="Model") +
    theme_bw() + theme(text = element_text("Calibri"),legend.position = "bottom") })
dev.off()


##Plot PPR posterior distribution
plotname <- ifelse(red==F,paste0("plots/mod_p_SMunif_",ctry,"_ppr_post_",joy,".png"),
                   paste0("plots/mod_p_SMunif_",ctry,"_ppr_post_",joy,"_red.png"))
png(file=plotname,width=9,height=10,units="cm",res=400,family="Calibri")
print({genunc$ppr %>% left_join(moddat) %>% left_join(pprobs,multiple="all") %>%
    left_join(data.frame(ind=2:5,ppr=c("0 -> 1","1 -> 2","2 -> 3","3+ -> 4+"))) %>%
    ggplot(aes(x=coh,color=factor(ppr),fill=factor(ppr))) +
    geom_ribbon(aes(ymin=lwr,ymax=upr,group=ppr),linetype=2,linewidth=0.3,alpha=0.3) +
    geom_line(aes(y=med,group=ppr)) +
    facet_grid(~paste0(ctryl,", ",joy),scales="free") +
    scale_x_continuous(breaks=seq(1900,2005,5),minor_breaks=1900:2100,expand=c(0.02,0.02)) +
    scale_y_continuous(breaks=seq(0,3,0.1)) +
    geom_point(aes(x=coh,y=obs,color=factor(ppr)),pch=21,fill="black",stroke=0.5,size=1,show.legend=F) +
    geom_vline(data=intdat %>% filter(joy2==joy),aes(xintercept=int)) +
    labs(x="Cohort",y="Posterior distribution",fill="PPR",color="PPR") +
    theme_bw() + theme(text = element_text("Calibri"),legend.position = "bottom") })
dev.off()


##Plot FSD posterior distribution
plotname <- ifelse(red==F,paste0("plots/mod_p_SMunif_",ctry,"_fsd_post_",joy,".png"),
                   paste0("plots/mod_p_SMunif_",ctry,"_fsd_post_",joy,"_red.png"))
png(file=plotname,width=8.5,height=10,units="cm",res=400,family="Calibri")
print({genunc$fsd %>% left_join(moddat) %>% left_join(fsdobs,multiple="all") %>%
    left_join(data.frame(ind=2:5,ppr=c("0","1","2","3+"))) %>%
    ggplot(aes(x=coh,color=factor(ppr),fill=factor(ppr))) +
    geom_ribbon(aes(ymin=lwr,ymax=upr,group=ppr),linetype=2,linewidth=0.3,alpha=0.3) +
    geom_line(aes(y=med,group=ppr)) +
    facet_grid(~paste0(ctryl,", ",joy),scales="free") +
    scale_x_continuous(breaks=seq(1900,2005,5),minor_breaks=1900:2100,expand=c(0.02,0.02)) +
    scale_y_continuous(breaks=seq(0,3,0.1)) +
    geom_point(aes(x=coh,y=obs,color=factor(ppr)),pch=21,fill="black",stroke=0.5,size=1,show.legend=F) +
    geom_vline(data=intdat %>% filter(joy2==joy),aes(xintercept=int)) +
    labs(x="Cohort",y="Posterior distribution",fill="No. of children",color="No. of children") +
    theme_bw() + theme(text = element_text("Calibri"),legend.position = "bottom") })
dev.off()
