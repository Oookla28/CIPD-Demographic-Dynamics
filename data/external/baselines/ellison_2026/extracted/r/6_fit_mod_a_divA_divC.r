###Fit model
##Load stan data
if (red==F) load(paste0("data/standata_",ctry,"_",joy,"_divA",divA,"_divC",divC,".RData"))
if (red==T) load(paste0("data/standata_",ctry,"_",joy,"_red_divA",divA,"_divC",divC,".RData"))


##Fit model
n_warm <- 250
n_keep <- 500
n_tot <- n_warm+n_keep
stanout <- stan(file="stan/mod_a.stan",
                data=standata,chains=1,iter=n_tot,warmup=n_warm)


##Save output
if (red==F) save(stanout,file=paste0("output/mod_a_",ctry,"_",joy,"_divA",divA,"_divC",divC,".RData"))
if (red==T) save(stanout,file=paste0("output/mod_a_",ctry,"_",joy,"_red_divA",divA,"_divC",divC,".RData"))


##Check convergence
plotname <- ifelse(red==F,paste0("plots/mod_a_",ctry,"_rhat_",joy,"_divA",divA,"_divC",divC,".png"),
                   paste0("plots/mod_a_",ctry,"_rhat_",joy,"_red_divA",divA,"_divC",divC,".png"))
png(file=plotname,width=15,height=10,units="cm",res=400,family="Calibri")
print({stan_rhat(stanout)})
dev.off()


##Inspect parameter samples (trace plots)
plotname <- ifelse(red==F,paste0("plots/mod_a_",ctry,"_trace1_",joy,"_divA",divA,"_divC",divC,".png"),
                   paste0("plots/mod_a_",ctry,"_trace1_",joy,"_red_divA",divA,"_divC",divC,".png"))
png(file=plotname,width=25,height=10,units="cm",res=400,family="Calibri")
print({stan_trace(stanout,pars=c("beta_1","beta_n",
                                 "sigmaA0r","sigmaC0r","phi0",
                                 paste0("beta0raw[",seq(4,standata$betadim-2*(standata$NdC+1),10),"]")),
                  inc_warmup = T,nrow=2)})
dev.off()
