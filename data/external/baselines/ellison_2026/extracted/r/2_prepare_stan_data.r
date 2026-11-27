###Prepare data for modelling
##Read in rates
load(paste0("data/",ctry,"_rates.RData"))


##Modelling setup
Na <- length(agerange) # number of ages (30)
Nc <- length(cohrange) # number of cohorts (40 or 35)
#Create data frame of all age-cohort combinations
dat <- expand.grid(Age=agerange,Cohort=cohrange)
Nf <- nrow(dat) # number of age-cohort combinations (30 x 40 or 35 = 1200 or 1050)
#Extract age-cohort exposures for full surface
ACexp <- ACexpmat$N[ACexpmat$Cohort %in% cohrange & ACexpmat$Age %in% agerange]

#Number of age-cohort combinations for each parity (Nf for ASFRs/parity 0, then
#Nc fewer for each successive parity, as we start from age 16 for parity 1, age 17
#for parity 2, and age 18 for parity 3+)
Nf1 <- Nf
Nf2 <- Nf-Nc
Nf3 <- Nf-2*Nc
Nf4 <- Nf-3*Nc

#Create indicators for age-cohort combinations omitted and included for parities 1, 2 and 3+
indA2 <- which(dat$Age<16)
indA2b <- setdiff(1:Nf,indA2)
indA3 <- which(dat$Age<17)
indA3b <- setdiff(1:Nf,indA3)
indA4 <- which(dat$Age<18)
indA4b <- setdiff(1:Nf,indA4)


##Prepare data
alldata$Year <- alldata$Age + alldata$Cohort #add year column
#All (ages 15-44)
data0 <- alldata[alldata$Year<=joy & alldata$Age %in% agerange & alldata$Cohort %in% cohrange & alldata$par==0,]
indp0 <- which(!is.na(left_join(data.frame(expand.grid(Age=agerange,Cohort=cohrange)),data0)$x))
N0 <- nrow(data0)
succ0 <- data0$n
tot0 <- data0$N

#Parity 0 (ages 15-44)
data1 <- alldata[alldata$Year<=joy & alldata$Age %in% agerange & alldata$Cohort %in% cohrange & alldata$par==1,]
indp1 <- which(!is.na(left_join(data.frame(expand.grid(Age=agerange,Cohort=cohrange)),data1)$par))
N1 <- nrow(data1)
succ1 <- data1$n
tot1 <- data1$N

#Parity 1 (ages 16-44)
data2 <- alldata[alldata$Year<=joy & alldata$Age %in% agerange[-c(1)] & alldata$Cohort %in% cohrange & alldata$par==2,]
indp2 <- which(!is.na(left_join(data.frame(expand.grid(Age=agerange,Cohort=cohrange)),data2)$par))
N2 <- nrow(data2)
succ2 <- data2$n
tot2 <- data2$N

#Parity 2 (ages 17-44)
data3 <- alldata[alldata$Year<=joy & alldata$Age %in% agerange[-c(1,2)] & alldata$Cohort %in% cohrange & alldata$par==3,]
indp3 <- which(!is.na(left_join(data.frame(expand.grid(Age=agerange,Cohort=cohrange)),data3)$par))
N3 <- nrow(data3)
succ3 <- data3$n
tot3 <- data3$N

#Parity 3+ (ages 18-44)
data4 <- alldata[alldata$Year<=joy & alldata$Age %in% agerange[-c(1,2,3)] & alldata$Cohort %in% cohrange & alldata$par==4,]
indp4 <- which(!is.na(left_join(data.frame(expand.grid(Age=agerange,Cohort=cohrange)),data4)$par))
N4 <- nrow(data4)
succ4 <- data4$n
tot4 <- data4$N


##Prepare basis functions and penalty matrices
NkA <- floor(Na/5) # number of knots for Age
NkC <- floor(Nc/5) # number of knots for cohort

#Construct cubic B-spline bases for Age
B_Af1 <- bbase(dat$Age, n.knots=NkA, deg=3)
NbA <- ncol(B_Af1) # basis dimension (NkA + 2)
NdA <- NbA-1 # number of first-order differences in age direction
B_Af2 <- bbase(dat$Age[dat$Age>agerange[1]], n.knots=NkA, deg=3)
B_Af3 <- bbase(dat$Age[dat$Age>agerange[2]], n.knots=NkA, deg=3)
B_Af4 <- bbase(dat$Age[dat$Age>agerange[3]], n.knots=NkA, deg=3)

#Construct cubic B-spline bases for cohort
B_Cf1 <- bbase(dat$Cohort, n.knots=NkC, deg=3)
NbC <- ncol(B_Cf1) # basis dimension (NkC + 2)
NdC <- NbC-1 # number of first-order differences in cohort direction
B_Cf2 <- bbase(dat$Cohort[dat$Age>agerange[1]], n.knots=NkC, deg=3)
B_Cf3 <- bbase(dat$Cohort[dat$Age>agerange[2]], n.knots=NkC, deg=3)
B_Cf4 <- bbase(dat$Cohort[dat$Age>agerange[3]], n.knots=NkC, deg=3)

#Compute number of first-order differences
NpA <- NdA*NbC # across age
NpC <- NdC*NbA # across cohort

#Compute number of tensor product basis functions
betadim <- NbA*NbC

#Construct tensor product B-spline bases 
B_ACf1 <- matrix(0,Nf,betadim)
B_ACf2 <- matrix(0,Nf-Nc,betadim)
B_ACf3 <- matrix(0,Nf-2*Nc,betadim)
B_ACf4 <- matrix(0,Nf-3*Nc,betadim)

k <- 1
for (i in 1:NbC) {
  for (j in 1:NbA) {
    B_ACf1[,k] <- B_Af1[,j]*B_Cf1[,i]
    B_ACf2[,k] <- B_Af2[,j]*B_Cf2[,i]
    B_ACf3[,k] <- B_Af3[,j]*B_Cf3[,i]
    B_ACf4[,k] <- B_Af4[,j]*B_Cf4[,i]
    k <- k+1 
  }
}

#Rename basis function matrices
Xf1 <- B_ACf1
Xf2 <- B_ACf2
Xf3 <- B_ACf3
Xf4 <- B_ACf4


##Construct difference matrices
#First-order difference matrix for age direction (1D)
DA <- diff(diag(NbA),1)
#Second-order difference matrix for age direction (1D)
DA2 <- diff(diff(diag(NbA),1),1)
#First-order difference matrix for cohort direction (1D)
DC <- diff(diag(NbC),1)
#Second-order difference matrix for cohort direction (1D)
DC2 <- diff(diff(diag(NbC),1),1)
#First order difference matrix for age direction (2D)
D1 <- kronecker(diag(NbC),DA)
#Second-order difference matrix for age direction (2D)
D12 <- kronecker(diag(NbC),DA2)
#First-order difference matrix for cohort direction (2D)
D2 <- kronecker(DC,diag(NbA))
#Second-order difference matrix for cohort direction (2D)
D22 <- kronecker(DC2,diag(NbA))

#Construct indicator matrix for computing CFRs from age-cohort surface of ASFRs
M <- matrix(0,Nc,Nf)
for (i in 1:Nc) M[i,(((i-1)*Na)+1):(i*Na)] <- 1


##Construct indicators
#Particular coefficients
ind1 <- seq(1,by=NbA,length=NbC) # first age row
indn <- seq(NbA,by=NbA,length=NbC) # last age row
indn1 <- seq(NbA-1,by=NbA,length=NbC) #penultimate age row
indr <- setdiff(1:betadim,c(ind1,indn)) # excluding first and last age rows
#First-order differences in cohort direction
ind21 <- seq(1,by=NbA,length=NbC-1) # first age row
ind2n <- (NbA-1)+ind21 # last age row
ind2n1 <- (NbA-2)+ind21 # penultimate age row
ind2r <- setdiff(1:NpC,c(ind21,ind2n)) # excluding first and last age rows
#First-order differences in age direction
ind11 <- seq(1,by=NbA-1,length=NbC) # first row of age differences
ind1n1 <- (NbA-3)+ind11 # penultimate row of age differences
ind1n <- (NbA-2)+ind11 # last row of age differences
ind1r <- setdiff(1:NpA,c(ind11,ind1n1,ind1n)) # excluding first, penultimate and last rows


##Standard deviation values for constraints
#Process HFD estimates of summary measures
source("r/1b_process_HFD_sum_data.r")
#PPRs
HFDdat5
pprsd <- signif(HFDdat5$sd,c(2,2,2,2))
#CFR
HFDdat6
cfrsd <- signif(HFDdat6$sd,2)


##Prior for PPR01 (Appendix F)
#Extract PPR01 values
ppr01 <- HFDdat3 %>% filter(name=="PPR01",!is.na(value)) %>% pull(value)

#Create function to calculate the sum of the absolute differences between the
#PPR01 density (for PPR01 values from the last local maximum onwards) and a lognormal CDF
#with parameters x[1] and x[2]
minfunc <- function(x) {
  dens <- density(ppr01)
  ind <<- length(dens$x)-(1:which(diff(rev(dens$y))<0)[1]-1)
  xval <- dens$x[ind]
  minx <<- min(xval)
  yval <- dens$y[ind]/dens$y[min(ind)]
  diff <- yval-(1-plnorm(xval,x[1],x[2]))
  sum(abs(diff))
}

#Identify the lognormal distribution parameters that minimise minfunc2
optout <- optim(c(log(0.95),0.25),minfunc)
ppr1_ln1b <- round(optout$par[1],3) #extract first parameter
ppr1_ln2b <- round(optout$par[2],3) #extract second parameter

#Plot density of PPR01 with lognormal CDF overlaid
png(file=paste0("plots/ppr01_prior_",joy,"_comp.png"),width=15,height=10,units="cm",res=400,family="Calibri")
func <- function(x,p1,p2) 1-plnorm(x,p1,p2)
dens <- density(ppr01)
df <- data.frame(x=dens$x,y=dens$y/max(dens$y[ind]))
print({ggplot() +
    geom_line(data=df,aes(x=x,y=y,color="Density (scaled)")) +
    scale_x_continuous(breaks=seq(0,1,0.02),limits=c(0.9,1)) +
    scale_y_continuous(breaks=seq(0,1,0.2),limits=c(0,1)) +
    geom_function(fun = func, args=list(p1=ppr1_ln1b,p2=ppr1_ln2b),aes(color="Lognormal CDF approximation")) +
    geom_vline(aes(xintercept=minx),linetype=2) +
    geom_hline(aes(yintercept=1),linetype=2) +
    scale_color_manual(values=c(1,2)) +
    labs(x="PPR 0 -> 1",y="Density",color="") +
    theme_bw() +
    theme(legend.position="bottom")})
dev.off()


##Stan setup
standata <- list(Nf1=Nf1,Nf2=Nf2,Nf3=Nf3,Nf4=Nf4,
                 Nc=Nc,Na=Na,NpA=NpA,NpC=NpC,NdC=NdC,
                 N0=N0,N1=N1,N2=N2,N3=N3,N4=N4,
                 ACexp=ACexp,betadim=betadim,M=M,
                 indA2=indA2,indA2b=indA2b,indA3=indA3,indA3b=indA3b,indA4=indA4,indA4b=indA4b,
                 ind1=ind1,indn=indn,indr=indr,ind1r=ind1r,ind2r=ind2r,
                 indp0=indp0,indp1=indp1,indp2=indp2,indp3=indp3,indp4=indp4,
                 succ0=succ0,succ1=succ1,succ2=succ2,succ3=succ3,succ4=succ4,
                 tot1=tot1,tot2=tot2,tot3=tot3,tot4=tot4,
                 Xf1=Xf1,Xf2=Xf2,Xf3=Xf3,Xf4=Xf4,
                 D1=D1,D2=D2,
                 pprsd=pprsd,cfrsd=cfrsd,
                 ppr1_ln1b=ppr1_ln1b,ppr1_ln2b=ppr1_ln2b)
if (red==F) save(standata,file=paste0("data/standata_",ctry,"_",joy,".RData"))
if (red==T) save(standata,file=paste0("data/standata_",ctry,"_",joy,"_red.RData"))
