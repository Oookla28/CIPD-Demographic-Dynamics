###Fit model
##Setup
n_warm <- 250
n_keep <- 500
n_tot <- n_warm+n_keep


##Load original model to provide initial values
if (red==F) load(file=paste0("output/mod_a_",ctry,"_",joy,".RData"))
if (red==T) load(file=paste0("output/mod_a_",ctry,"_",joy,"_red.RData"))

inits <- list()
inits[[1]] <- list()
inits[[1]]$beta_1 <- extract(stanout,"beta_1",permuted=F)[,1,] %>% median()
inits[[1]]$beta_n <- extract(stanout,"beta_n",permuted=F)[,1,] %>% median()
inits[[1]]$beta0raw <- extract(stanout,"beta0raw",permuted=F)[,1,] %>% apply(2,median)
inits[[1]]$sigmaA0r <- extract(stanout,"sigmaA0r",permuted=F)[,1,] %>% median()
inits[[1]]$sigmaC0r <- extract(stanout,"sigmaC0r",permuted=F)[,1,] %>% median()
inits[[1]]$sigmaCFR <- mean(standata$cfr_lu)
inits[[1]]$phi0 <- extract(stanout,"phi0",permuted=F)[,1,] %>% median()


##Fit model
stanout <- stan(file="stan/mod_a_SMunif.stan",
                data=standata,chains=1,iter=n_tot,warmup=n_warm,
                init = inits)


##Save output
if (red==F) save(stanout,file=paste0("output/mod_a_SMunif_",ctry,"_",joy,".RData"))
if (red==T) save(stanout,file=paste0("output/mod_a_SMunif_",ctry,"_",joy,"_red.RData"))


##Check convergence
plotname <- ifelse(red==F,paste0("plots/mod_a_SMunif_",ctry,"_rhat_",joy,".png"),
                   paste0("plots/mod_a_SMunif_",ctry,"_rhat_",joy,"_red.png"))
png(file=plotname,width=15,height=10,units="cm",res=400,family="Calibri")
print({stan_rhat(stanout)})
dev.off()


##Inspect parameter samples (trace plots)
plotname <- ifelse(red==F,paste0("plots/mod_a_SMunif_",ctry,"_trace1_",joy,".png"),
                   paste0("plots/mod_a_SMunif_",ctry,"_trace1_",joy,"_red.png"))
png(file=plotname,width=25,height=10,units="cm",res=400,family="Calibri")
print({stan_trace(stanout,pars=c("beta_1","beta_n",
                                 "sigmaA0r","sigmaC0r","sigmaCFR","phi0",
                                 paste0("beta0raw[",seq(4,54,10),"]")),
                  inc_warmup = T,nrow=2)})
dev.off()
