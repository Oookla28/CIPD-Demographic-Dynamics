functions {
  real lognorm_pseudo_log(vector x, real ln1, real ln2) {
    vector[num_elements(x)] lp1;
    real lp;
    real tmp;
    for (i in 1:num_elements(x)) {
      tmp = 1-lognormal_cdf(x[i],ln1,ln2)+0.00001;
      lp1[i] = log(tmp);
    }
    lp = sum(lp1);
    return lp;
  }
}

data {
  int Nf1;
  int Nf2;
  int Nf3;
  int Nf4;
  int Nc;
  int Nc_f1;
  int Na;
  int Na_obs[Nc];
  int NpA;
  int NpC;
  int NdC;
  int N0;
  int N1;
  int N2;
  int N3;
  int N4;  
  vector[Nf1] ACexp;
  int betadim;
  matrix[Nc,Nf1] M;
  int indA2[Nc];
  int indA2b[Nf1-Nc];
  int indA3[2*Nc];
  int indA3b[Nf1-2*Nc];
  int indA4[3*Nc];
  int indA4b[Nf1-3*Nc];
  int ind1[NdC+1];
  int indn[NdC+1];
  int indr[betadim-2*(NdC+1)];
  int ind1r[NpA-3*NdC-3];
  int ind2r[NpC-2*NdC];
  int indp0[N0];
  int indp1[N1];
  int indp2[N2];
  int indp3[N3];
  int indp4[N4];
  int succ0[N0];
  int succ1[N1];
  int succ2[N2];
  int succ3[N3];
  int succ4[N4];
  vector[N1] tot1;
  vector[N2] tot2;
  vector[N3] tot3;
  vector[N4] tot4;
  matrix[N0,4] tot_pac;
  matrix[Nf1,betadim] Xf1;
  matrix[Nf2,betadim] Xf2;
  matrix[Nf3,betadim] Xf3;
  matrix[Nf4,betadim] Xf4;
  matrix[NpA,betadim] D1;
  matrix[NpC,betadim] D2;
  real pprsd[4];
  real cfrsd;
  real ppr1_ln1b;
  real ppr1_ln2b;
}

parameters {
  vector[4] beta_1;
  vector[4] beta_n;
  vector[betadim-2*(NdC+1)] beta1raw;
  vector[betadim-2*(NdC+1)] beta2raw;
  vector[betadim-2*(NdC+1)] beta3raw;
  vector[betadim-2*(NdC+1)] beta4raw;
  real<lower=0> sigmaA1r;
  real<lower=0> sigmaC1r;
  real<lower=0> sigmaA2r;
  real<lower=0> sigmaC2r;
  real<lower=0> sigmaA3r;
  real<lower=0> sigmaC3r;
  real<lower=0> sigmaA4r;
  real<lower=0> sigmaC4r;
  real<lower=0> phi0;
  real<lower=0> phi1;
  real<lower=0> phi2;
  real<lower=0> phi3;
  real<lower=0> phi4;
}

transformed parameters {
  vector[Nf1] Y_hatf1;
  vector[Nf1] Y_hatf2;
  vector[Nf1] Y_hatf3;
  vector[Nf1] Y_hatf4;
  vector[betadim] beta1;
  vector[betadim] beta2;
  vector[betadim] beta3;
  vector[betadim] beta4;
  vector[Nf1] pf1;
  vector[Nf1] pf2;
  vector[Nf1] pf3;
  vector[Nf1] pf4;
  vector[Nc] cfr;
  vector[Nc] cfr1;
  vector[Nc] cfr2;
  vector[Nc] cfr3;
  vector[Nc] cfr4;
  vector[Nc] ppr1;
  vector[Nc] ppr2;
  vector[Nc] ppr3;
  vector[Nc] ppr4;
  vector[Nf1] btot;
  vector[Nf1] asfr;
  vector[Nf1] asfr1;
  vector[Nf1] asfr2;
  vector[Nf1] asfr3;
  vector[Nf1] asfr4;
  matrix[Nf1,4] births;
  matrix[Nf1,4] expos;
  row_vector[4] nxtn;
  vector[NpA-3*NdC-3] D1rbeta1;
  vector[NpC-2*NdC] D2rbeta1;
  vector[NpA-3*NdC-3] D1rbeta2;
  vector[NpC-2*NdC] D2rbeta2;
  vector[NpA-3*NdC-3] D1rbeta3;
  vector[NpC-2*NdC] D2rbeta3;
  vector[NpA-3*NdC-3] D1rbeta4;
  vector[NpC-2*NdC] D2rbeta4;
  beta1[indr] = beta1raw;
  beta2[indr] = beta2raw;
  beta3[indr] = beta3raw;
  beta4[indr] = beta4raw;
  beta1[ind1] = rep_vector(beta_1[1],NdC+1);
  beta2[ind1] = rep_vector(beta_1[2],NdC+1);
  beta3[ind1] = rep_vector(beta_1[3],NdC+1);
  beta4[ind1] = rep_vector(beta_1[4],NdC+1);
  beta1[indn] = rep_vector(beta_n[1],NdC+1);
  beta2[indn] = rep_vector(beta_n[2],NdC+1);
  beta3[indn] = rep_vector(beta_n[3],NdC+1);
  beta4[indn] = rep_vector(beta_n[4],NdC+1);
  D2rbeta1 = D2[ind2r,]*beta1;
  D1rbeta1 = D1[ind1r,]*beta1;
  D2rbeta2 = D2[ind2r,]*beta2;
  D1rbeta2 = D1[ind1r,]*beta2;
  D2rbeta3 = D2[ind2r,]*beta3;
  D1rbeta3 = D1[ind1r,]*beta3;
  D2rbeta4 = D2[ind2r,]*beta4;
  D1rbeta4 = D1[ind1r,]*beta4;
  Y_hatf1 = Xf1*beta1;
  Y_hatf2[indA2] = rep_vector(0,Nc);
  Y_hatf2[indA2b] = Xf2*beta2;
  Y_hatf3[indA3] = rep_vector(0,2*Nc);
  Y_hatf3[indA3b] = Xf3*beta3;
  Y_hatf4[indA4] = rep_vector(0,3*Nc);
  Y_hatf4[indA4b] = Xf4*beta4;
  pf1 = exp(Y_hatf1);
  pf2 = exp(Y_hatf2);
  pf2[indA2] = rep_vector(0,Nc);
  pf3 = exp(Y_hatf3);
  pf3[indA3] = rep_vector(0,2*Nc);
  pf4 = exp(Y_hatf4);
  pf4[indA4] = rep_vector(0,3*Nc);
  //Exposures and births for observed age-cohort combinations
  for (i in 1:N0) {
    expos[indp0[i],] = tot_pac[i,];
    births[indp0[i],1] = expos[indp0[i],1]*pf1[indp0[i]];
    births[indp0[i],2] = expos[indp0[i],2]*pf2[indp0[i]];
    births[indp0[i],3] = expos[indp0[i],3]*pf3[indp0[i]];
    births[indp0[i],4] = expos[indp0[i],4]*pf4[indp0[i]];
  }
  //Exposures and births for future age-cohort combinations
  for (i in Nc_f1:Nc) {
    for (j in Na_obs[i]:Na) {
      if (j > Na_obs[i]) {
        births[(i-1)*Na+j,1] = expos[(i-1)*Na+j,1]*pf1[(i-1)*Na+j];
        births[(i-1)*Na+j,2] = expos[(i-1)*Na+j,2]*pf2[(i-1)*Na+j];
        births[(i-1)*Na+j,3] = expos[(i-1)*Na+j,3]*pf3[(i-1)*Na+j];
        births[(i-1)*Na+j,4] = expos[(i-1)*Na+j,4]*pf4[(i-1)*Na+j];
      }
      if (j < Na) {
        nxtn[1] = expos[(i-1)*Na+j,1]-births[(i-1)*Na+j,1];
        if (nxtn[1] < 0) {
          nxtn[1] = 1;
        }
        nxtn[2] = expos[(i-1)*Na+j,2]-births[(i-1)*Na+j,2]+births[(i-1)*Na+j,1];
        if (nxtn[2] < 0) {
          nxtn[2] = 1;
        }
        nxtn[3] = expos[(i-1)*Na+j,3]-births[(i-1)*Na+j,3]+births[(i-1)*Na+j,2];
        if (nxtn[3] < 0) {
          nxtn[3] = 1;
        }
        nxtn[4] = expos[(i-1)*Na+j,4]+births[(i-1)*Na+j,3];
        if (nxtn[4] < 0) {
          nxtn[4] = 1;
        }
        expos[(i-1)*Na+j+1,] = ACexp[(i-1)*Na+j+1]*nxtn/sum(nxtn);
      }
    }
  }
  for (i in 1:Nf1) {
    btot[i] = sum(births[i,]);
    asfr[i] = btot[i]/ACexp[i];
    asfr1[i] = births[i,1]/ACexp[i];
    asfr2[i] = births[i,2]/ACexp[i];
    asfr3[i] = births[i,3]/ACexp[i];
    asfr4[i] = births[i,4]/ACexp[i];
  }
  cfr = M * asfr;
  cfr1 = M * asfr1;
  cfr2 = M * asfr2;
  cfr3 = M * asfr3;
  cfr4 = M * asfr4;
  ppr1 = cfr1;
  for (i in 1:Nc) {
    ppr2[i] = cfr2[i]/cfr1[i];
    ppr3[i] = cfr3[i]/cfr2[i];
    ppr4[i] = cfr4[i]/(cfr3[i]+cfr4[i]);
  }
}

model {
  D2rbeta1 ~ normal(0,sigmaC1r);
  D1rbeta1 ~ normal(0,sigmaA1r);
  D2rbeta2 ~ normal(0,sigmaC2r);
  D1rbeta2 ~ normal(0,sigmaA2r);
  D2rbeta3 ~ normal(0,sigmaC3r);
  D1rbeta3 ~ normal(0,sigmaA3r);
  D2rbeta4 ~ normal(0,sigmaC4r);
  D1rbeta4 ~ normal(0,sigmaA4r);
  sigmaC1r ~ student_t(7,0,0.01);
  sigmaA1r ~ student_t(7,0,0.01);
  sigmaC2r ~ student_t(7,0,0.01);
  sigmaA2r ~ student_t(7,0,0.01);
  sigmaC3r ~ student_t(7,0,0.01);
  sigmaA3r ~ student_t(7,0,0.01);
  sigmaC4r ~ student_t(7,0,0.01);
  sigmaA4r ~ student_t(7,0,0.01);
  phi0 ~ normal(0,1000);
  phi1 ~ normal(0,5);
  phi2 ~ normal(0,5);
  phi3 ~ normal(0,5);
  phi4 ~ normal(0,5);
  for (n in 2:Nc) {
    ppr1[n] ~ normal(ppr1[n-1],pprsd[1]);
    ppr2[n] ~ normal(ppr2[n-1],pprsd[2]);
    ppr3[n] ~ normal(ppr3[n-1],pprsd[3]);
    ppr4[n] ~ normal(ppr4[n-1],pprsd[4]);
    cfr[n] ~ normal(cfr[n-1],cfrsd);
  }
  ppr1 ~ lognorm_pseudo(ppr1_ln1b,ppr1_ln2b);
  for (n in 1:N0) {
   succ0[n] ~ neg_binomial_2(btot[indp0[n]],phi0);
  }
  for (n in 1:N1) {
    succ1[n] ~ neg_binomial_2(tot1[n]*pf1[indp1[n]],phi1);
  }
  for (n in 1:N2) {
    succ2[n] ~ neg_binomial_2(tot2[n]*pf2[indp2[n]],phi2);
  }
  for (n in 1:N3) {
    succ3[n] ~ neg_binomial_2(tot3[n]*pf3[indp3[n]],phi3);
  }
  for (n in 1:N4) {
    succ4[n] ~ neg_binomial_2(tot4[n]*pf4[indp4[n]],phi4);
  }  
}
