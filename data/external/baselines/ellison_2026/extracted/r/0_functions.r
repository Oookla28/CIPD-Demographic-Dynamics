###Important functions
##Functions taken from Miller and Neocleous (2018) Flexible Regression Lecture Notes 
tpower <- function(x, t, p)
  (x - t) ^ p * (x > t)

bbase <- function(x, xl = min(x), xr = max(x), n.knots = 10, deg = 3) {
  # Construct B-spline basis (based on a function written by Paul Eilers)
  nseg <- n.knots-1
  dx <- (xr - xl) / nseg
  knots <- seq(xl - deg * dx, xr + deg * dx, len = n.knots + 2*deg )
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
  B <- (-1) ^ (deg + 1) * P %*% t(D)
  B
}


##Function to generate uncertainty for an aggregate or parity-specific model
#Arguments:
#files    - vector of file locations of stan outputs of models
#agerange - range of age variable
#joyvec   - corresponding vector of jump-off years
#modvec   - corresponding vector of model types ("Aggregate" or "Parity-specific")
#ctryvec  - corresponding vector of countries
#redvec   - corresponding vector of indicators (FALSE if the full 40 cohorts are modelled,
#           TRUE if only the reduced 35 cohorts are modelled)
#regexp   - corresponding vector of indicators (FALSE if all parity-specific exposures are
#           approximated using the cumulation-based approach, TRUE if register-based exposures
#           are used for the observed age-cohort combinations)
gen_unc <- function(files,agerange,
                    joyvec=rep(joy,length(files)),
                    modvec=rep(mods,each=length(ctrys)),
                    ctryvec=rep(ctrys,2),
                    redvec=rep(red,length(files)),
                    regexp=rep(F,length(files))) {
  for (f in 1:length(files)) {
    cat(paste0("File ",f," of ",length(files),"\n"))
    
    #Setup
    if (redvec[f]==F) load(paste0("data/standata_",ctryvec[f],"_",joyvec[f],".RData"))
    if (redvec[f]==T) load(paste0("data/standata_",ctryvec[f],"_",joyvec[f],"_red.RData"))
    attach(standata, warn.conflicts = F)
    load(paste0("output/",files[f]))
    if (redvec[f]==F) cohrange <- (joyvec[f]-max(agerange)-10):(joyvec[f]-min(agerange))
    if (redvec[f]==T) cohrange <- (joyvec[f]-max(agerange)-5):(joyvec[f]-min(agerange))
    Nc <- length(cohrange)
    Na <- length(agerange)
    Nf <- Na*Nc
    Ns <- as.numeric(unlist(stanout@stan_args)["iter"])-as.numeric(unlist(stanout@stan_args)["warmup"])
    
    #Aggregate level
    phi0 <- extract(stanout,"phi0",permuted=F)[,1,]
    phi0 <- as.matrix(phi0)[,rep(1,Nf)]
    asfr <- extract(stanout,"asfr",permuted=F)[,1,]
    
    #Generate births
    set.seed(1)
    ACexp <- standata$ACexp
    births_agg <- matrix(
      rnbinom(prod(Nf,Ns),
              mu=asfr*matrix(rep(ACexp,Ns),nrow=Ns,byrow=T),size=phi0),nrow=Ns,byrow=F)
    
    #Summarise posterior distribution of ASFRs
    asfrdat <- data.frame(mod=f,expand.grid(age=agerange,coh=cohrange),med=NA,lwr=NA,upr=NA)
    asfrdat$med <- apply(asfr,2,median)
    asfrdat$lwr <- apply(births_agg,2,quantile,p=0.025)/ACexp
    asfrdat$upr <- apply(births_agg,2,quantile,p=0.975)/ACexp
    asfrdat$mod <- factor(asfrdat$mod)
    
    #Summarise posterior distribution of CFR
    M <- matrix(0,Nc,Nf)
    for (i in 1:Nc) M[i,(((i-1)*Na)+1):(i*Na)] <- 1
    
    cfrdat <- data.frame(mod=f,coh=cohrange,ind=1,med=NA,lwr=NA,upr=NA)
    cfrdat$med <- apply((births_agg/matrix(rep(ACexp,Ns),nrow=Ns,byrow=T)) %*% t(M),2,median)
    cfrdat$lwr <- apply((births_agg/matrix(rep(ACexp,Ns),nrow=Ns,byrow=T)) %*% t(M),2,quantile,p=0.025)
    cfrdat$upr <- apply((births_agg/matrix(rep(ACexp,Ns),nrow=Ns,byrow=T)) %*% t(M),2,quantile,p=0.975)
    cfrdat$mod <- factor(cfrdat$mod)
    
    #Parity-specific level
    if (modvec[f]=="Parity-specific") {
      phi <- extract(stanout,c("phi1","phi2","phi3","phi4"),permuted=F)[,1,]
      phi <- phi[,c(rep(1,Nf),rep(2,Nf),rep(3,Nf),rep(4,Nf))]
      pf <- extract(stanout,c("pf1","pf2","pf3","pf4"),permuted=F)[,1,]
      btot <- matrix(extract(stanout,"btot",permuted=F)[,1,],nrow=Ns,byrow=F)
      pf[,c(1*Nf+standata$indA2,2*Nf+standata$indA3,3*Nf+standata$indA4)] <- 0
      
      if (regexp[f] == F) {
        #Generate order-specific births and compute denominators
        births_par <- expos_par <- matrix(0,Ns,4*Nf)
        nxtn <- numeric(4)
        for (k in 1:Ns) {
          if(k %% 100 == 0) cat(100*k/Ns,"%\n",sep="")
          for (i in 1:Nc) {
            for (j in 1:Na) {
              if (j == 1) {
                expos_par[k,0*Nf+(i-1)*Na+j] = ACexp[(i-1)*Na+j];
                expos_par[k,1*Nf+(i-1)*Na+j] = 0;
                expos_par[k,2*Nf+(i-1)*Na+j] = 0;
                expos_par[k,3*Nf+(i-1)*Na+j] = 0;
              }
              if (expos_par[k,0*Nf+(i-1)*Na+j]<0) expos_par[k,0*Nf+(i-1)*Na+j] <- 0
              if (expos_par[k,1*Nf+(i-1)*Na+j]<0) expos_par[k,1*Nf+(i-1)*Na+j] <- 0
              if (expos_par[k,2*Nf+(i-1)*Na+j]<0) expos_par[k,2*Nf+(i-1)*Na+j] <- 0
              if (expos_par[k,3*Nf+(i-1)*Na+j]<0) expos_par[k,3*Nf+(i-1)*Na+j] <- 0
              births_par[k,0*Nf+(i-1)*Na+j] = rnbinom(1,mu=expos_par[k,0*Nf+(i-1)*Na+j]*pf[k,0*Nf+(i-1)*Na+j],size=phi[k,0*Nf+(i-1)*Na+j]);
              births_par[k,1*Nf+(i-1)*Na+j] = rnbinom(1,mu=expos_par[k,1*Nf+(i-1)*Na+j]*pf[k,1*Nf+(i-1)*Na+j],size=phi[k,1*Nf+(i-1)*Na+j]);
              births_par[k,2*Nf+(i-1)*Na+j] = rnbinom(1,mu=expos_par[k,2*Nf+(i-1)*Na+j]*pf[k,2*Nf+(i-1)*Na+j],size=phi[k,2*Nf+(i-1)*Na+j]);      
              births_par[k,3*Nf+(i-1)*Na+j] = rnbinom(1,mu=expos_par[k,3*Nf+(i-1)*Na+j]*pf[k,3*Nf+(i-1)*Na+j],size=phi[k,3*Nf+(i-1)*Na+j]);
              if (j < Na) {
                nxtn[1] = expos_par[k,0*Nf+(i-1)*Na+j]-births_par[k,0*Nf+(i-1)*Na+j];
                nxtn[2] = expos_par[k,1*Nf+(i-1)*Na+j]-births_par[k,1*Nf+(i-1)*Na+j]+births_par[k,0*Nf+(i-1)*Na+j];
                nxtn[3] = expos_par[k,2*Nf+(i-1)*Na+j]-births_par[k,2*Nf+(i-1)*Na+j]+births_par[k,1*Nf+(i-1)*Na+j];
                nxtn[4] = expos_par[k,3*Nf+(i-1)*Na+j]+births_par[k,2*Nf+(i-1)*Na+j];
                expos_par[k,c(0*Nf,1*Nf,2*Nf,3*Nf)+(i-1)*Na+j+1] = ACexp[(i-1)*Na+j+1]*nxtn/sum(nxtn);
              }
            }
          }
        }
      }
      
      if (regexp[f] == T) {
        #Load standata
        if (redvec[f]==F) load(paste0("data/standata_",ctryvec[f],"_",joyvec[f],"_pac.RData"))
        if (redvec[f]==T) load(paste0("data/standata_",ctryvec[f],"_",joyvec[f],"_pac_red.RData"))
        #Generate order-specific births and compute denominators for future age-cohort combinations only
        births_par <- expos_par <- matrix(0,Ns,4*Nf)
        nxtn <- numeric(4)
        for (k in 1:Ns) {
          if (k %% 100 == 0) cat(100*k/Ns,"%\n",sep="")
          #Observed age-cohort combinations
          expos_par[k,0*Nf+standata$indp0] = standata$tot_pac[,1]
          expos_par[k,1*Nf+standata$indp0] = standata$tot_pac[,2]
          expos_par[k,2*Nf+standata$indp0] = standata$tot_pac[,3]
          expos_par[k,3*Nf+standata$indp0] = standata$tot_pac[,4]
          for (i in 1:standata$N0) {
            births_par[k,0*Nf+standata$indp0[i]] =  rnbinom(1,mu=expos_par[k,0*Nf+standata$indp0[i]]*pf[k,0*Nf+standata$indp0[i]],size=phi[k,0*Nf+standata$indp0[i]]);
            births_par[k,1*Nf+standata$indp0[i]] =  rnbinom(1,mu=expos_par[k,1*Nf+standata$indp0[i]]*pf[k,1*Nf+standata$indp0[i]],size=phi[k,1*Nf+standata$indp0[i]]);
            births_par[k,2*Nf+standata$indp0[i]] =  rnbinom(1,mu=expos_par[k,2*Nf+standata$indp0[i]]*pf[k,2*Nf+standata$indp0[i]],size=phi[k,2*Nf+standata$indp0[i]]);
            births_par[k,3*Nf+standata$indp0[i]] =  rnbinom(1,mu=expos_par[k,3*Nf+standata$indp0[i]]*pf[k,3*Nf+standata$indp0[i]],size=phi[k,3*Nf+standata$indp0[i]]);
          }
          #Future age-cohort combinations
          for (i in standata$Nc_f1:Nc) {
            for (j in standata$Na_obs[i]:Na) {
              if (expos_par[k,0*Nf+(i-1)*Na+j]<0) expos_par[k,0*Nf+(i-1)*Na+j] <- 0
              if (expos_par[k,1*Nf+(i-1)*Na+j]<0) expos_par[k,1*Nf+(i-1)*Na+j] <- 0
              if (expos_par[k,2*Nf+(i-1)*Na+j]<0) expos_par[k,2*Nf+(i-1)*Na+j] <- 0
              if (expos_par[k,3*Nf+(i-1)*Na+j]<0) expos_par[k,3*Nf+(i-1)*Na+j] <- 0
              if (j > standata$Na_obs[i]) {
                births_par[k,0*Nf+(i-1)*Na+j] = rnbinom(1,mu=expos_par[k,0*Nf+(i-1)*Na+j]*pf[k,0*Nf+(i-1)*Na+j],size=phi[k,0*Nf+(i-1)*Na+j]);
                births_par[k,1*Nf+(i-1)*Na+j] = rnbinom(1,mu=expos_par[k,1*Nf+(i-1)*Na+j]*pf[k,1*Nf+(i-1)*Na+j],size=phi[k,1*Nf+(i-1)*Na+j]);
                births_par[k,2*Nf+(i-1)*Na+j] = rnbinom(1,mu=expos_par[k,2*Nf+(i-1)*Na+j]*pf[k,2*Nf+(i-1)*Na+j],size=phi[k,2*Nf+(i-1)*Na+j]);      
                births_par[k,3*Nf+(i-1)*Na+j] = rnbinom(1,mu=expos_par[k,3*Nf+(i-1)*Na+j]*pf[k,3*Nf+(i-1)*Na+j],size=phi[k,3*Nf+(i-1)*Na+j]);
              }              
              if (j < Na) {
                nxtn[1] = expos_par[k,0*Nf+(i-1)*Na+j]-births_par[k,0*Nf+(i-1)*Na+j];
                nxtn[2] = expos_par[k,1*Nf+(i-1)*Na+j]-births_par[k,1*Nf+(i-1)*Na+j]+births_par[k,0*Nf+(i-1)*Na+j];
                nxtn[3] = expos_par[k,2*Nf+(i-1)*Na+j]-births_par[k,2*Nf+(i-1)*Na+j]+births_par[k,1*Nf+(i-1)*Na+j];
                nxtn[4] = expos_par[k,3*Nf+(i-1)*Na+j]+births_par[k,2*Nf+(i-1)*Na+j];
                expos_par[k,c(0*Nf,1*Nf,2*Nf,3*Nf)+(i-1)*Na+j+1] = ACexp[(i-1)*Na+j+1]*nxtn/sum(nxtn);
              }
            }
          }
        }
      }
      
      rates_par <- births_par/expos_par
      expos_par_t <- expos_par[,1:Nf]+expos_par[,Nf+1:Nf]+expos_par[,2*Nf+1:Nf]+expos_par[,3*Nf+1:Nf]
      
      #Summarise posterior distribution of conditional ASFRs
      casfrdat <- data.frame(mod=f,expand.grid(age=agerange,coh=cohrange,par=1:4),med=NA,lwr=NA,upr=NA)
      casfrdat$med <- apply(pf,2,median,na.rm=T)
      casfrdat$lwr <- apply(rates_par,2,quantile,p=0.025,na.rm=T)
      casfrdat$upr <- apply(rates_par,2,quantile,p=0.975,na.rm=T)
      casfrdat$mod <- factor(casfrdat$mod)
      
      #Compute unconditional ASFRs
      uasfr1 <- births_par[,1:Nf]/expos_par_t
      uasfr2 <- births_par[,Nf+1:Nf]/expos_par_t
      uasfr3 <- births_par[,2*Nf+1:Nf]/expos_par_t
      uasfr4 <- births_par[,3*Nf+1:Nf]/expos_par_t
      
      #Compute summary measures
      cfr1 <- uasfr1 %*% t(M)
      cfr2 <- uasfr2 %*% t(M)
      cfr3 <- uasfr3 %*% t(M)
      cfr4 <- uasfr4 %*% t(M)
      
      ppr1 <- cfr1
      ppr2 <- cfr2/cfr1
      ppr3 <- cfr3/cfr2
      ppr4 <- cfr4/(cfr3+cfr4)
      pprall <- cbind(ppr1,ppr2,ppr3,ppr4)
      
      fsd1 <- 1-cfr1
      fsd2 <- cfr1-cfr2
      fsd3 <- cfr2-cfr3
      fsd4 <- cfr3
      fsdall <- cbind(fsd1,fsd2,fsd3,fsd4)
      
      #Summarise posterior distribution of summary measures
      pprdat <- data.frame(mod=f,expand.grid(coh=cohrange,ind=2:5),med=NA,lwr=NA,upr=NA)
      pprdat$med <- apply(pprall,2,median)
      pprdat$lwr <- apply(pprall,2,quantile,p=0.025)
      pprdat$upr <- apply(pprall,2,quantile,p=0.975)
      pprdat$mod <- factor(pprdat$mod)
      
      cfrdat2 <- data.frame(mod=f,coh=cohrange,ind=1,med=NA,lwr=NA,upr=NA)
      cfrdat2$med <- apply(((births_par[,1:Nf]+births_par[,Nf+1:Nf]+births_par[,2*Nf+1:Nf]+births_par[,3*Nf+1:Nf])/matrix(rep(ACexp,Ns),nrow=Ns,byrow=T)) %*% t(M),2,median)
      cfrdat2$lwr <- apply(((births_par[,1:Nf]+births_par[,Nf+1:Nf]+births_par[,2*Nf+1:Nf]+births_par[,3*Nf+1:Nf])/matrix(rep(ACexp,Ns),nrow=Ns,byrow=T)) %*% t(M),2,quantile,p=0.025)
      cfrdat2$upr <- apply(((births_par[,1:Nf]+births_par[,Nf+1:Nf]+births_par[,2*Nf+1:Nf]+births_par[,3*Nf+1:Nf])/matrix(rep(ACexp,Ns),nrow=Ns,byrow=T)) %*% t(M),2,quantile,p=0.975)
      cfrdat2$mod <- factor(cfrdat2$mod)
      
      fsddat <- data.frame(mod=f,expand.grid(coh=cohrange,ind=2:5),med=NA,lwr=NA,upr=NA)
      fsddat$med <- apply(fsdall,2,median)
      fsddat$lwr <- apply(fsdall,2,quantile,p=0.025)
      fsddat$upr <- apply(fsdall,2,quantile,p=0.975)
      fsddat$mod <- factor(fsddat$mod)
    }
    
    #Extract output
    if (f==1) {
      out <- list()
      out$asfr  <- asfrdat
      out$cfr   <- cfrdat
      if (modvec[f]=="Parity-specific") out$casfr <- casfrdat
      if (modvec[f]=="Parity-specific") out$ppr   <- pprdat
      if (modvec[f]=="Parity-specific") out$cfr2  <- cfrdat2
      if (modvec[f]=="Parity-specific") out$fsd   <- fsddat
    }
    if (f>1) {
      out$asfr  <- rbind(out$asfr,asfrdat)
      out$cfr   <- rbind(out$cfr,cfrdat)
      if (modvec[f]=="Parity-specific") out$casfr <- rbind(out$casfr,casfrdat)
      if (modvec[f]=="Parity-specific") out$ppr   <- rbind(out$ppr,pprdat)
      if (modvec[f]=="Parity-specific") out$cfr2  <- rbind(out$cfr2,cfrdat2)
      if (modvec[f]=="Parity-specific") out$fsd   <- rbind(out$fsd,fsddat)
    }
  }
  return(out)
}


##Function to generate detailed aggregate level summary statistics for an aggregate or parity-specific model
#Arguments:
#files    - vector of file locations of stan outputs of models
#agerange - range of age variable
#joyvec   - corresponding vector of jump-off years
#modvec   - corresponding vector of model types ("Aggregate" or "Parity-specific")
#ctryvec  - corresponding vector of countries
#redvec   - corresponding vector of indicators (FALSE if the full 40 cohorts are modelled,
#           TRUE if only the reduced 35 cohorts are modelled)
gen_unc2 <- function(files,agerange,
                     joyvec=rep(joy,length(files)),
                     modvec=rep(mods,each=length(ctrys)),
                     ctryvec=rep(ctrys,2),
                     redvec=rep(red,length(files))) {
  for (f in 1:length(files)) {
    cat(paste0("File ",f," of ",length(files),"\n"))
    
    #Setup
    if (redvec[f]==F) load(paste0("data/standata_",ctryvec[f],"_",joyvec[f],".RData"))
    if (redvec[f]==T) load(paste0("data/standata_",ctryvec[f],"_",joyvec[f],"_red.RData"))
    attach(standata, warn.conflicts = F)
    load(paste0("output/",files[f]))
    if (redvec[f]==F) cohrange <- (joyvec[f]-max(agerange)-10):(joyvec[f]-min(agerange))
    if (redvec[f]==T) cohrange <- (joyvec[f]-max(agerange)-5):(joyvec[f]-min(agerange))
    Nc <- length(cohrange)
    Na <- length(agerange)
    Nf <- Na*Nc
    Ns <- as.numeric(unlist(stanout@stan_args)["iter"])-as.numeric(unlist(stanout@stan_args)["warmup"])
    
    #Process data
    obsasfr <- alldata2 %>% filter(ctry==ctryvec[f],par==0,Cohort %in% cohrange) %>%
      select(Cohort,Age,Year,x) %>% rename("age"=Age,"coh"=Cohort,"yr"=Year)
    coh.comp <- as.numeric(names(which(table(obsasfr$coh)==Na)))
    obscfr <- obsasfr %>% filter(coh %in% coh.comp, yr<=joyvec[f]) %>%
      group_by(coh) %>% summarise(jcfr=sum(x))
    obscfr2 <- obsasfr %>% filter(coh %in% coh.comp) %>%
      group_by(coh) %>% summarise(ocfr=sum(x))
    
    #Aggregate level
    phi0 <- extract(stanout,"phi0",permuted=F)[,1,]
    phi0 <- as.matrix(phi0)[,rep(1,Nf)]
    asfr <- extract(stanout,"asfr",permuted=F)[,1,]
    
    #Generate births
    set.seed(1)
    ACexp <- standata$ACexp
    births_agg <- matrix(
      rnbinom(prod(Nf,Ns),
              mu=asfr*matrix(rep(ACexp,Ns),nrow=Ns,byrow=T),size=phi0),nrow=Ns,byrow=F)
    
    #Summarise posterior distribution of ASFRs
    asfrdat <- data.frame(mod=f,expand.grid(age=agerange,coh=cohrange))
    asfrdat$yr <- asfrdat$age+asfrdat$coh
    asfrdat$med <- apply(births_agg,2,median)/ACexp
    asfrdat$mean <- apply(births_agg,2,mean)/ACexp
    asfrdat$sd <- apply(births_agg,2,sd)/ACexp
    asfrdat$lwr95 <- apply(births_agg,2,quantile,p=0.025)/ACexp
    asfrdat$upr95 <- apply(births_agg,2,quantile,p=0.975)/ACexp
    asfrdat$lwr90 <- apply(births_agg,2,quantile,p=0.05)/ACexp
    asfrdat$upr90 <- apply(births_agg,2,quantile,p=0.95)/ACexp
    asfrdat$lwr50 <- apply(births_agg,2,quantile,p=0.25)/ACexp
    asfrdat$upr50 <- apply(births_agg,2,quantile,p=0.75)/ACexp
    asfrdat$mod <- factor(asfrdat$mod)
    asfrdat %<>% left_join(obsasfr)
    
    #Prepare for CRPS calculation
    ind <- which(asfrdat$yr>joyvec[f] & !is.na(asfrdat$x))
    tmp2 <- t(t(births_agg)/ACexp)
    
    #Calculate CRPS (using the empirical CDF-based approximation given in equation (3) of Jordan et al. (2019), DOI: 10.18637/jss.v090.i12)
    asfrdat$crps <- NA
    for (i in ind) {
      asfrdat$crps[i] <- (2/(Ns^2))*sum((tmp2[,i] - asfrdat$x[i])*(Ns*as.numeric(rep(asfrdat$x[i], Ns) < tmp2[,i]) - (1:Ns) + 0.5))
    }
    
    #Calculate LogS (under a Gaussian approximation)
    asfrdat$logs <- NA
    asfrdat$logs[ind] <- -log(dnorm(asfrdat$x[ind],asfrdat$mean[ind],asfrdat$sd[ind]))
    
    #Summarise posterior distribution of CFR
    M <- matrix(0,Nc,Nf)
    for (i in 1:Nc) M[i,(((i-1)*Na)+1):(i*Na)] <- 1
    M2 <- M
    M2[,which(asfrdat$age+asfrdat$coh<=joyvec[f])] <- 0
    
    #Summing across all ages
    cfrdat <- data.frame(mod=f,coh=cohrange)
    cfrdat$med <- apply(t(t(births_agg)/ACexp) %*% t(M),2,median)
    cfrdat$mean <- apply(t(t(births_agg)/ACexp) %*% t(M),2,mean)
    cfrdat$sd <- apply(t(t(births_agg)/ACexp) %*% t(M),2,sd)
    cfrdat$lwr95 <- apply(t(t(births_agg)/ACexp) %*% t(M),2,quantile,p=0.025)
    cfrdat$upr95 <- apply(t(t(births_agg)/ACexp) %*% t(M),2,quantile,p=0.975)
    cfrdat$lwr90 <- apply(t(t(births_agg)/ACexp) %*% t(M),2,quantile,p=0.05)
    cfrdat$upr90 <- apply(t(t(births_agg)/ACexp) %*% t(M),2,quantile,p=0.95)
    cfrdat$lwr50 <- apply(t(t(births_agg)/ACexp) %*% t(M),2,quantile,p=0.25)
    cfrdat$upr50 <- apply(t(t(births_agg)/ACexp) %*% t(M),2,quantile,p=0.75)
    cfrdat$mod <- factor(cfrdat$mod)
    
    #Summing across the ages lived by the cohort after the jump-off year
    jcfrdat <- data.frame(mod=f,coh=cohrange)
    jcfrdat$med <- apply(t(t(births_agg)/ACexp) %*% t(M2),2,median)
    jcfrdat$mean <- apply(t(t(births_agg)/ACexp) %*% t(M2),2,mean)
    jcfrdat$sd <- apply(t(t(births_agg)/ACexp) %*% t(M2),2,sd)
    jcfrdat$lwr95 <- apply(t(t(births_agg)/ACexp) %*% t(M2),2,quantile,p=0.025)
    jcfrdat$upr95 <- apply(t(t(births_agg)/ACexp) %*% t(M2),2,quantile,p=0.975)
    jcfrdat$lwr90 <- apply(t(t(births_agg)/ACexp) %*% t(M2),2,quantile,p=0.05)
    jcfrdat$upr90 <- apply(t(t(births_agg)/ACexp) %*% t(M2),2,quantile,p=0.95)
    jcfrdat$lwr50 <- apply(t(t(births_agg)/ACexp) %*% t(M2),2,quantile,p=0.25)
    jcfrdat$upr50 <- apply(t(t(births_agg)/ACexp) %*% t(M2),2,quantile,p=0.75)
    jcfrdat$mod <- factor(jcfrdat$mod)
    jcfrdat %<>% left_join(obscfr %>% left_join(obscfr2))
    
    #Prepare for CRPS calculation
    jcfrdat %<>% mutate(med2=jcfr+med,mean2=jcfr+mean,
                        lwr95_2=jcfr+lwr95,upr95_2=jcfr+upr95,
                        lwr90_2=jcfr+lwr90,upr90_2=jcfr+upr90,
                        lwr50_2=jcfr+lwr50,upr50_2=jcfr+upr50)
    ind <- which(cohrange+max(agerange)>joyvec[f] & !is.na(jcfrdat$ocfr))
    tmp2 <- t(t(births_agg)/ACexp) %*% t(M2)
    
    #Calculate CRPS (using the empirical CDF-based approximation given in equation (3) of Jordan et al. (2019), DOI: 10.18637/jss.v090.i12)
    jcfrdat$crps <- NA
    for (i in ind) {
      jcfrdat$crps[i] <- (2/(Ns^2))*sum((tmp2[,i] + jcfrdat$jcfr[i] - jcfrdat$ocfr[i])*(Ns*as.numeric(rep(jcfrdat$ocfr[i], Ns) < tmp2[,i] + jcfrdat$jcfr[i]) - (1:Ns) + 0.5))
    }
    
    #Calculate LogS (under a Gaussian approximation)
    jcfrdat$logs <- NA
    jcfrdat$logs[ind] <- -log(dnorm(jcfrdat$ocfr[ind],jcfrdat$mean2[ind],jcfrdat$sd[ind]))
    
    #Extract output
    if (f==1) {
      out <- list()
      out$asfr  <- asfrdat
      out$cfr   <- cfrdat
      out$jcfr  <- jcfrdat
    }
    
    if (f>1) {
      out$asfr  <- rbind(out$asfr,asfrdat)
      out$cfr   <- rbind(out$cfr,cfrdat)
      out$jcfr  <- rbind(out$jcfr,jcfrdat)
    }
  }
  return(out)
}


##Function to calculate CRPS for a Gaussian predictive distribution
#(see p367 of Gneiting and Raftery (2007), DOI: 10.1198/016214506000001437)
#Arguments:
#x  - true value
#mu - mean of predictive distribution
#sd - standard deviation of predictive distribution
CRPS_normal_func <- function(x, mu, sigma) {
  -sigma*((1/sqrt(pi)) - 2*dnorm((x - mu)/sigma) - ((x-mu)/sigma)*(2*pnorm((x - mu)/sigma) - 1))
}
