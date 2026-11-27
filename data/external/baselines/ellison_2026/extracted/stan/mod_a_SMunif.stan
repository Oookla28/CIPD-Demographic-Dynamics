data {
  int Nf1;
  int Nc;
  int NpA;
  int NpC;
  int NdC;
  int N0;
  vector[Nf1] ACexp;
  int betadim;
  matrix[Nc,Nf1] M;
  int ind1[NdC+1];
  int indn[NdC+1];
  int indr[betadim-2*(NdC+1)];
  int ind1r[NpA-3*NdC-3];
  int ind2r[NpC-2*NdC];
  int indp0[N0];
  int succ0[N0];
  matrix[Nf1,betadim] Xf1;
  matrix[NpA,betadim] D1;
  matrix[NpC,betadim] D2;
  real cfrsd;
  vector[2] cfr_lu;
}

parameters {
  real beta_1;
  real beta_n;
  vector[betadim-2*(NdC+1)] beta0raw;
  real<lower=0> sigmaA0r;
  real<lower=0> sigmaC0r;
  real<lower=cfr_lu[1],upper=cfr_lu[2]> sigmaCFR;
  real<lower=0> phi0;
}

transformed parameters {
  vector[Nf1] Y_hatf0;
  vector[betadim] beta0;
  vector[Nf1] asfr;
  vector[Nc] cfr;
  vector[NpA-3*NdC-3] D1rbeta0;
  vector[NpC-2*NdC] D2rbeta0;
  beta0[indr] = beta0raw;
  beta0[ind1] = rep_vector(beta_1,NdC+1);
  beta0[indn] = rep_vector(beta_n,NdC+1);
  D2rbeta0 = D2[ind2r,]*beta0;
  D1rbeta0 = D1[ind1r,]*beta0;
  Y_hatf0 = Xf1*beta0;
  asfr = exp(Y_hatf0);
  cfr = M * asfr;
}

model {
  D2rbeta0 ~ normal(0,sigmaC0r);
  D1rbeta0 ~ normal(0,sigmaA0r);
  sigmaC0r ~ student_t(7,0,0.01);
  sigmaA0r ~ student_t(7,0,0.01);
  phi0 ~ normal(0,1000);
  for (n in 2:Nc) {
    cfr[n] ~ normal(cfr[n-1],sigmaCFR);
  }
  sigmaCFR ~ uniform(cfr_lu[1],cfr_lu[2]);
  for (n in 1:N0) {
   succ0[n] ~ neg_binomial_2(ACexp[indp0[n]]*asfr[indp0[n]],phi0);
  }
}
