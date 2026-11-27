###Compare summary measure SD posterior distributions
load("output/resSMunif.RData")


##Prepare results objects
Niter <- 500 #number of iterations
resdat <- matrix(0,Niter,6*3*length(joys))
ctry_vec <- numeric()
joy_vec <- numeric()


##Extract model output
k <- 1
for (i in 1:length(joys)) { #for each JOY
  #Set JOY
  joy <- joys[i]
  print(joy)
  
  #Identify countries
  ctrys_unif <- resSMunif[[i]]$ctry$Code
  print(ctrys_unif)
  
  for (ctry in ctrys_unif) { #for each country
    #Setup
    print(ctry)
    print(k)
    joy_vec[k] <- joy
    ctry_vec[k] <- ctry
    joy_sel <- joy_sel_lst[[ctry]]
    red_sel <- red_sel_lst[[ctry]]
    red <- red_sel[which(joy_sel==joy)]
    
    #Load aggregate model output
    if (red==F) load(paste0("output/mod_a_SMunif_",ctry,"_",joy,".RData")); stanouta <- stanout
    if (red==T) load(paste0("output/mod_a_SMunif_",ctry,"_",joy,"_red.RData")); stanouta <- stanout
    
    #Load parity-specific model output
    if (red==F) load(paste0("output/mod_p_SMunif_",ctry,"_",joy,".RData")); stanoutp <- stanout
    if (red==T) load(paste0("output/mod_p_SMunif_",ctry,"_",joy,"_red.RData")); stanoutp <- stanout
    
    #Extract parameter estimates
    resdat[,(k-1)*6+1:6] <- cbind(extract(stanouta, pars="sigmaCFR",permuted=F)[,1,],
                                   extract(stanoutp, pars=c("sigmaCFR",paste0("sigmaPPR",1:4)),
                                           permuted=F)[,1,])

    k <- k+1
  }
}


##Plot densities
#Prepare plotting object
resdat2 <- resdat %>% cbind(iter=1:Niter) %>% as_tibble() %>% 
  pivot_longer(-iter,names_to="par") %>%
  mutate(par=as.numeric(substr(par,2,10))) %>%
  mutate(par2=rep(paste0("SM = ",c("CFR (A)","CFR (P-S)",paste0("PPR ",c(0:2,"3+")," -> ",c(1:3,"4+")))),3*length(joys)*Niter)) %>%
  mutate(ctry=ctry_vec[ceiling(par/6)],
         joy=joy_vec[ceiling(par/6)])

#Extract SD boundaries and average SD for plotting
bnd_avg_dat <- lapply(resSMunif,function(x) x$bounds[c(1:2,1:10),c("sd","sdt")]) %>% bind_rows() %>%
  mutate(joy=rep(joys,each=12),
         par2=rep(rep(paste0("SM = ",c("CFR (A)","CFR (P-S)",paste0("PPR ",c(0:2,"3+")," -> ",c(1:3,"4+")))),each=2),length(joys)))

#Save plot
png(file=paste0("plots/mod_ap_SMunif_sd_dens_comp.png"),width=31,height=15,units="cm",res=400,family="Calibri")
print({resdat2 %>% ggplot(aes(x=value,color=ctry)) +
  geom_density() +
  geom_vline(data=bnd_avg_dat,aes(xintercept=sd),lty=1) +
  geom_vline(data=bnd_avg_dat,aes(xintercept=sdt),lty=2) +
  facet_wrap(~paste0("JOY = ",joy,", ",par2),ncol=6,scales="free") +
  labs(x=expression(sigma),y="Density",color="Country") +
  theme_bw() +
  guides(color=guide_legend(nrow=1)) +
  theme(legend.position = "bottom")})
dev.off()
