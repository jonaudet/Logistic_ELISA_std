data{
int N_unkn;           //Total number of unknown values
int N_unkn_grp;       //Total number of unknown samples (should be N_unkn / 4)
int uID[N_unkn];      //Sample number for each unknown value
vector[N_unkn] Unknown; //Optical Density of each unknown
vector[N_unkn] ser_dilutions; //The serial dilutions from the start dilution for each unknown;
  real log_mu_Std;
  real<lower = 0> log_sigma_std;
}
transformed data{
  vector[N_unkn] log_ser_dilutions;
  vector[N_unkn] abs_log_ser_dilutions;
  int std_loc;

  log_ser_dilutions <- log(ser_dilutions);
  std_loc <- max(uID);
  for(i in 1:N_unkn)
    abs_log_ser_dilutions[i] <- fabs(log_ser_dilutions[i]);
}
parameters{
  real<lower = 0> sigma;


  real mu_Span;
  real mu_Bottom;
  real mu_log_Inflec;
  real<lower = 0> mu_Slope;


  vector[N_unkn] log_x_raw; // log of each unknown's initial concentration
  vector[N_unkn_grp - 1] log_theta;  // log of each unknown's predicted concentration
  real<lower = 0> sigma_x;
  real pred_std_raw;
}
model{
  vector[N_unkn] unkn_cOD;
  vector[N_unkn] log_x;
  vector[N_unkn] log_Undil;
  real log_pred_std;

  sigma ~ normal(0, 1);
  sigma_x ~ normal(0, 1);
  mu_Bottom ~ normal(0.05, 0.01);
  mu_Span ~ normal(3.5, 0.1);
  mu_log_Inflec ~ normal(0, 2);
  mu_Slope ~ normal(1, 0.5);

  //Multilevel unknown estimation
  log_theta ~ normal(5, 10);
  log_x_raw ~ normal(0, 1);
  pred_std_raw ~ normal(0, 1);

  log_pred_std <- log_mu_Std + pred_std_raw * log_sigma_std;

  for(i in 1:N_unkn){
    if(uID[i] == std_loc){
      log_Undil[i] <- log_pred_std;
    } else {
      log_Undil[i] <- log_theta[uID[i]];
    }
  }

  log_x <- (log_ser_dilutions + log_Undil) + (sigma_x * abs_log_ser_dilutions) .* log_x_raw;

  for(i in 1:N_unkn){
    unkn_cOD[i] <- mu_Bottom + mu_Span * inv_logit((log_x[i] - mu_log_Inflec) * mu_Slope);
  }

  Unknown ~ normal(unkn_cOD, sigma);
}
generated quantities{
  vector[N_unkn_grp - 1] theta;
  vector[N_unkn] x;

{
  vector[N_unkn] log_Undil;
  real log_pred_std;

  log_pred_std <- log_mu_Std + pred_std_raw * log_sigma_std;

  for(i in 1:N_unkn){
    if(uID[i] == std_loc){
      log_Undil[i] <- log_pred_std;
    } else {
      log_Undil[i] <- log_theta[uID[i]];
    }
  }

  x <- exp((log_ser_dilutions + log_Undil) + (sigma_x * abs_log_ser_dilutions) .* log_x_raw);
}
  theta <- exp(log_theta);
}
