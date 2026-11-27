###Fit model
##Setup
n_warm <- 250
n_keep <- 500
n_tot <- n_warm+n_keep


##Load original model to provide initial values
if (red==F) load(file=paste0("output/mod_p_",ctry,"_",joy,".RData"))
if (red==T) load(file=paste0("output/mod_p_",ctry,"_",joy,"_red.RData"))

inits <- list()
inits[[1]] <- list()
inits[[1]]$beta_1 <- extract(stanout,"beta_1",permuted=F)[,1,] %>% apply(2,median)
inits[[1]]$beta_n <- extract(stanout,"beta_n",permuted=F)[,1,] %>% apply(2,median)
inits[[1]]$beta1raw <- extract(stanout,"beta1raw",permuted=F)[,1,] %>% apply(2,median)
inits[[1]]$beta2raw <- extract(stanout,"beta2raw",permuted=F)[,1,] %>% apply(2,median)
inits[[1]]$beta3raw <- extract(stanout,"beta3raw",permuted=F)[,1,] %>% apply(2,median)
inits[[1]]$beta4raw <- extract(stanout,"beta4raw",permuted=F)[,1,] %>% apply(2,median)
inits[[1]]$sigmaA1r <- extract(stanout,"sigmaA1r",permuted=F)[,1,] %>% median()
inits[[1]]$sigmaC1r <- extract(stanout,"sigmaC1r",permuted=F)[,1,] %>% median()
inits[[1]]$sigmaA2r <- extract(stanout,"sigmaA2r",permuted=F)[,1,] %>% median()
inits[[1]]$sigmaC2r <- extract(stanout,"sigmaC2r",permuted=F)[,1,] %>% median()
inits[[1]]$sigmaA3r <- extract(stanout,"sigmaA3r",permuted=F)[,1,] %>% median()
inits[[1]]$sigmaC3r <- extract(stanout,"sigmaC3r",permuted=F)[,1,] %>% median()
inits[[1]]$sigmaA4r <- extract(stanout,"sigmaA4r",permuted=F)[,1,] %>% median()
inits[[1]]$sigmaC4r <- extract(stanout,"sigmaC4r",permuted=F)[,1,] %>% median()
inits[[1]]$sigmaCFR <- mean(standata$cfr_lu)
inits[[1]]$sigmaPPR <- apply(matrix(standata$ppr_lu,ncol=2,byrow=T),1,mean)
inits[[1]]$phi0 <- extract(stanout,"phi0",permuted=F)[,1,] %>% median()
inits[[1]]$phi1 <- extract(stanout,"phi1",permuted=F)[,1,] %>% median()
inits[[1]]$phi2 <- extract(stanout,"phi2",permuted=F)[,1,] %>% median()
inits[[1]]$phi3 <- extract(stanout,"phi3",permuted=F)[,1,] %>% median()
inits[[1]]$phi4 <- extract(stanout,"phi4",permuted=F)[,1,] %>% median()


##Fit model
stanout <- stan(file="stan/mod_p_SMunif.stan",
                data=standata,chains=1,iter=n_tot,warmup=n_warm,
                init = inits)


##Save output
if (red==F) save(stanout,file=paste0("output/mod_p_SMunif_",ctry,"_",joy,".RData"))
if (red==T) save(stanout,file=paste0("output/mod_p_SMunif_",ctry,"_",joy,"_red.RData"))


##Check convergence
plotname <- ifelse(red==F,paste0("plots/mod_p_SMunif_",ctry,"_rhat_",joy,".png"),
                   paste0("plots/mod_p_SMunif_",ctry,"_rhat_",joy,"_red.png"))
png(file=plotname,width=15,height=10,units="cm",res=400,family="Calibri")
print({stan_rhat(stanout)})
dev.off()


##Inspect parameter samples (trace plots)
plotname <- ifelse(red==F,paste0("plots/mod_p_SMunif_",ctry,"_trace1_",joy,".png"),
                   paste0("plots/mod_p_SMunif_",ctry,"_trace1_",joy,"_red.png"))
png(file=plotname,width=25,height=10,units="cm",res=400,family="Calibri")
print({stan_trace(stanout,pars=c("beta_1","beta_n",
                          paste0("sigmaA",1:4,"r"),
                          paste0("sigmaC",1:4,"r"),
                          "sigmaCFR","sigmaPPR1","sigmaPPR2","sigmaPPR3","sigmaPPR4",
                          paste0("phi",0:4)),inc_warmup = T,nrow=3)})
dev.off()

plotname <- ifelse(red==F,paste0("plots/mod_p_SMunif_",ctry,"_trace2_",joy,".png"),
                   paste0("plots/mod_p_SMunif_",ctry,"_trace2_",joy,"_red.png"))
png(file=plotname,width=25,height=10,units="cm",res=400,family="Calibri")
print({stan_trace(stanout,pars=paste0("beta",rep(1:4,each=6),
                                      "raw[",rep(seq(4,54,10),4),"]"),
                  inc_warmup = T,nrow=3)})
dev.off()
