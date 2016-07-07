data{
int N_unkn;           //Total number of unknown values
int N_unkn_grp;       //Total number of unknown samples (should be N_unkn / 4)
int uID[N_unkn];      //Sample number for each unknown value
vector[N_unkn] OD; //Optical Density to be modeled
vector[N_unkn] ser_dilutions; //The serial dilutions from the start dilution for each unknown;
}
transformed data{
  vector[N_unkn] log_ser_dilutions;

  log_ser_dilutions <- log(ser_dilutions);
}
parameters{
  real<lower = 0> sigma;


  real mu_Span;
  real mu_Bottom;
  real mu_log_Inflec;
  real<lower = 0> sigma_log_Inflec;
  vector[N_unkn_grp] log_Inflec_raw;
  real<upper = 0> mu_Slope;
}
transformed parameters{
  vector[N_unkn_grp] log_Inflec;

  log_Inflec <- mu_log_Inflec + sigma_log_Inflec * log_Inflec_raw;
}
model{
  vector[N_unkn] unkn_cOD;
  vector[N_unkn] log_x;
  vector[N_unkn] log_Undil;
  real pred_std;

  sigma ~ normal(0, 1);
  mu_Bottom ~ normal(0.05, 0.01);
  mu_Span ~ normal(3.5, 0.1);
  mu_log_Inflec ~ normal(0, 2);
  sigma_log_Inflec ~ normal(0, 1);
  log_Inflec_raw ~ normal(0, 1);
  mu_Slope ~ normal(-1, 0.5);

  for(i in 1:N_unkn){
    unkn_cOD[i] <- mu_Bottom + mu_Span * inv_logit((log_ser_dilutions[i] - log_Inflec[uID[i]]) * mu_Slope);
  }

  OD ~ normal(unkn_cOD, sigma);
}
