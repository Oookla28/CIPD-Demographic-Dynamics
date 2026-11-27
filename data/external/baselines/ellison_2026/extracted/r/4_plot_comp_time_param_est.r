###Plot computation times and parameter estimates
##Prepare data frames to store outputs
timedat <- sigmadat <- phidat <- data.frame()
sigmatmp <- data.frame(mod=rep(c("a","p"),c(2,8)),
                       cov=c("a","c",rep("a",4),rep("c",4)),
                       par=c(rep(0,2),rep(1:4,2)))
phitmp <- data.frame(mod=rep(c("a","p"),c(1,5)),
                     par=c(0,0:4))


##Extract information from fitted models
for (ctry in ctrys_fit) {
  print(ctry)
  joy_sel <- joy_sel_lst[[ctry]]
  red_sel <- red_sel_lst[[ctry]]
  for (i in 1:length(joy_sel)) {
    joy <- joy_sel[i]; print(joy)
    red <- red_sel[i]
    if (red==F) load(paste0("output/mod_a_",ctry,"_",joy,".RData"))
    if (red==T) load(paste0("output/mod_a_",ctry,"_",joy,"_red.RData"))
    stanouta <- stanout
    if (red==F) load(paste0("output/mod_p_",ctry,"_",joy,".RData"))
    if (red==T) load(paste0("output/mod_p_",ctry,"_",joy,"_red.RData"))
    stanoutp <- stanout
    sigma <- cbind(extract(stanouta,c("sigmaA0r","sigmaC0r"),permuted=F)[,1,],
                   extract(stanoutp,c(paste0("sigmaA",1:4,"r"),
                                      paste0("sigmaC",1:4,"r")),permuted=F)[,1,])
    phi <- cbind(extract(stanouta,c("phi0"),permuted=F)[,1,],
                 extract(stanoutp,c(paste0("phi",0:4)),permuted=F)[,1,])
    timedat <- rbind(timedat,data.frame(ctry=ctry,joy=joy,red=red,
                                        mod=c("a","p"),
                                        total=c(get_elapsed_time(stanouta) %>% sum(),
                                                get_elapsed_time(stanoutp) %>% sum()),
                                        warmup=c(get_elapsed_time(stanouta)[1],
                                                 get_elapsed_time(stanoutp)[1]),
                                        sample=c(get_elapsed_time(stanouta)[2],
                                                 get_elapsed_time(stanoutp)[2])))
    sigmadat <- rbind(sigmadat,data.frame(ctry=ctry,joy=joy,red=red,sigmatmp,
                                          med=apply(sigma,2,median),
                                          lwr=apply(sigma,2,quantile,p=0.025),
                                          upr=apply(sigma,2,quantile,p=0.975)))
    phidat <- rbind(phidat,data.frame(ctry=ctry,joy=joy,red=red,phitmp,
                                      med=apply(phi,2,median),
                                      lwr=apply(phi,2,quantile,p=0.025),
                                      upr=apply(phi,2,quantile,p=0.975)))
  }
}

#Compute average times for both model types
timedat %>% group_by(mod) %>% summarise(avg=mean(total)/60)


##Plot computation times
png(file=paste0("plots/mod_ap_comp_time.png"),width=12.5,height=10,units="cm",res=400,family="Calibri")
print({timedat %>% mutate(mod=ifelse(mod=="a","Aggregate","Parity-specific")) %>%
    ggplot(aes(y=ctry,x=total/60,color=factor(joy),shape=mod)) +
    geom_point() +
    labs(y="Country",x="Total time (minutes)",color="JOY",shape="Model") +
    scale_shape_manual(values=c(1,19)) +
    theme_bw() +
    theme(legend.position="right")})
dev.off()


##Plot parameter estimates
#Sigma
supp.labp <- c("Aggregate",paste0("Parity ",c(0,1,2,"3+")))
names(supp.labp) <- 0:4
supp.labc <- c("Age direction","Cohort direction"); names(supp.labc) <- c("a","c")

png(file="plots/mod_ap_sigma_par_est.png",width=25,height=15,units="cm",res=400,family="Calibri")
print({sigmadat %>% ggplot(aes(x=ctry,y=med, color=as.factor(joy))) +
  facet_grid(cov~par,labeller=labeller(par=supp.labp,cov=supp.labc)) +
  geom_errorbar(aes(ymin=lwr,ymax=upr)) +
  geom_point() +
  theme_bw() +
  labs(x="Country",y="Sigma",color="JOY") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position = "bottom")})
dev.off()

#Phi
supp.labp <- c("Aggregate",paste0("Parity ",c(0,1,2,"3+")))
names(supp.labp) <- 0:4

png(file="plots/mod_ap_phi_par_est.png",width=30,height=10,units="cm",res=400,family="Calibri")
print({phidat %>% 
    mutate(mod=ifelse(mod=="a","Aggregate","Parity-specific")) %>%
    ggplot(aes(x=ctry,y=med, color=as.factor(joy), shape=as.factor(mod))) +
    facet_wrap(~par,labeller=labeller(par=supp.labp),scales="free_y",nrow=1) +
    geom_errorbar(aes(ymin=lwr,ymax=upr)) +
    geom_point() +
    #coord_cartesian(ylim=c(0,1000)) +
    theme_bw() +
    scale_shape_manual(values=c(1,19)) +
    labs(x="Country",y="Phi",color="JOY",shape="Model") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))})
dev.off()
