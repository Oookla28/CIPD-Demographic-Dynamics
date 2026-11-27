###Fit model
##Load stan data
if (red==F) load(paste0("data/standata_",ctry,"_",joy,".RData"))
if (red==T) load(paste0("data/standata_",ctry,"_",joy,"_red.RData"))


##Fit model
n_warm <- 250
n_keep <- 500
n_tot <- n_warm+n_keep
stanout <- stan(file="stan/mod_p.stan",
                data=standata,chains=1,iter=n_tot,warmup=n_warm)


##Save output
if (red==F) save(stanout,file=paste0("output/mod_p_",ctry,"_",joy,".RData"))
if (red==T) save(stanout,file=paste0("output/mod_p_",ctry,"_",joy,"_red.RData"))


##Check convergence
plotname <- ifelse(red==F,paste0("plots/mod_p_",ctry,"_rhat_",joy,".png"),
                   paste0("plots/mod_p_",ctry,"_rhat_",joy,"_red.png"))
png(file=plotname,width=15,height=10,units="cm",res=400,family="Calibri")
print({stan_rhat(stanout)})
dev.off()


##Inspect parameter samples (trace plots)
plotname <- ifelse(red==F,paste0("plots/mod_p_",ctry,"_trace1_",joy,".png"),
                   paste0("plots/mod_p_",ctry,"_trace1_",joy,"_red.png"))
png(file=plotname,width=25,height=10,units="cm",res=400,family="Calibri")
print({stan_trace(stanout,pars=c("beta_1","beta_n",
                          paste0("sigmaA",1:4,"r"),
                          paste0("sigmaC",1:4,"r"),
                          paste0("phi",0:4)),inc_warmup = T,nrow=3)})
dev.off()

plotname <- ifelse(red==F,paste0("plots/mod_p_",ctry,"_trace2_",joy,".png"),
                   paste0("plots/mod_p_",ctry,"_trace2_",joy,"_red.png"))
png(file=plotname,width=25,height=10,units="cm",res=400,family="Calibri")
print({stan_trace(stanout,pars=paste0("beta",rep(1:4,each=6),
                                      "raw[",rep(seq(4,54,10),4),"]"),
                  inc_warmup = T,nrow=3)})
dev.off()
