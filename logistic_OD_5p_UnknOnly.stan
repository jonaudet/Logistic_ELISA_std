data{
int N_unkn;           //Total number of unknown values
int N_unkn_grp;       //Total number of unknown samples (should be N_unkn / 4)
int uID[N_unkn];      //Sample number for each unknown value
vector[N_unkn] Unknown; //Optical Density of each unknown
vector[N_unkn] ser_dilutions; //The serial dilutions from the start dilution for each unknown;
int pID_unkn[N_unkn]; //plate # for each unknown (to match it to the proper standard)
}
transformed data{
vector[N_unkn] log_ser_dilutions;

log_ser_dilutions <- log(ser_dilutions);
}
parameters{
  real<lower = 0> sigma;       //StDev of average measurement, as per Gelman et al

  
  real mu_Span;
  real mu_Bottom;
  real mu_log_Inflec;
  real<lower = 0> mu_Slope;


  vector[N_unkn] log_x_init; // log of each unknown's initial concentration
  vector[N_unkn_grp] log_theta;  // log of each unknown's predicted concentration
  vector<lower = 0>[N_unkn_grp] sigma_init;
  // real<lower = 0> avg_unkn_sigma;
}
transformed parameters{
// vector[N] sigma_indiv;         //The standard deviation for each point
vector[N_unkn] unkn_cOD;
// vector[N_unkn] unkn_sigma;
real mu_Inflec;
// vector[N_plates] Asyms;
vector[N_unkn] x_init;
// vector[N_unkn] x_i_unkn;


mu_Inflec <- exp(mu_log_Inflec);

//Unknown expected OD

x_init <- exp(log_x_init);

for(i in 1:N_unkn){
  // unkn_cOD[i] <- Bottoms[pID_unkn[i]] + (Spans[pID_unkn[i]] / (1 + (Inflecs[pID_unkn[i]] / x_i_unkn[i]) ^ Slopes[pID_unkn[i]]) ^ Asyms[pID_unkn[i]]);
  unkn_cOD[i] <- mu_Bottom + (mu_Span / (1 + (mu_Inflec / x_init[i]) ^ mu_Slope));
  // unkn_sigma[i] <- ((unkn_cOD[i] / 2) ^ (2 * alpha_unkn)) * avg_unkn_sigma;
}
}
model{
  sigma ~ cauchy(0, 1);
  mu_Bottom ~ normal(0.05, 0.01);
  mu_Span ~ normal(3.5, 0.1);
  mu_log_Inflec ~ uniform(-5, 10);
  mu_Slope ~ normal(1, 0.5);

  //Multilevel unknown estimation
  log_theta ~ uniform(-10, 15);
  log_x_init ~ normal(log_ser_dilutions + log_theta[uID], sigma_init[uID]);
  sigma_init ~ exponential(4);

  Unknown ~ normal(unkn_cOD, sigma);
}
generated quantities{
  vector[N_unkn_grp] theta;
  // real mu_x_i_unkn;
//   real sigma_Bottom;
//   real sigma_Span;
//   real sigma_Inflec;
//   real sigma_Slope;
//   real sigma_Asym;
//   real mu_Bottom;
//   real mu_Span;
//   real mu_Inflec;
//   real mu_Slope;
//   real mu_Asym;
//   real cv_Bottom;
//   real cv_Span;
//   real cv_Inflec;
//   real cv_Slope;
//   real cv_Asym;
// 
//   sigma_Bottom <- sd(Bottoms);
//   sigma_Span <- sd(Spans);
//   sigma_Inflec <- sd(Inflecs);
//   sigma_Slope <- sd(Slopes);
//   sigma_Asym <- sd(Asyms);
// 
//   mu_Bottom <- mean(Bottoms);
//   mu_Span <- mean(Spans);
//   mu_Inflec <- mean(Inflecs);
//   mu_Slope <- mean(Slopes);
//   mu_Asym <- mean(Asyms);
// 
//   cv_Bottom <- sigma_Bottom / mu_Bottom * 100;
//   cv_Span <- sigma_Span / mu_Span * 100;
//   cv_Inflec <- sigma_Inflec / mu_Inflec * 100;
//   cv_Slope <- sigma_Slope / mu_Slope * 100;
//   cv_Asym <- sigma_Asym / mu_Asym * 100;
  // mu_x_i_unkn <- mean(x_i_unkn);
  
  theta <- exp(log_theta);
}