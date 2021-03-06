data{
  int N;      //Total number of wells measured
  int N_grp;  //Total number of independent samples (1 sample has multiple dilutions which have multiple replicates)
  int N_grp_dil;  //Total number of dilutions and groups (so that N_grp_dil X # of replicates == N)
  int N_bot; //Number of values at Conc = 0

  int<lower = 1, upper = N_grp> uID[N_grp_dil]; //Sample ID for each dilution
  int<lower = 1, upper = N_grp_dil> dil_ID[N];  //Dilution ID for each well

  vector[N] meas_Signal;    //Optical density measured for each well
  vector[N_grp_dil] dilution;  //Dilution factor for each dilution

  real mu_Std;       //Concentration of the standard
  real<lower = 0> sigma_std;    //Uncertainty on the concentration of the standard

  vector[N_bot] zeroes;

  real inflec_mu;   //
}
transformed data{
  vector[N_grp_dil] log_dilution;
  vector[N_grp_dil] abs_log_dilution;
  int std_loc;

  log_dilution = log(dilution);
  for(i in 1:N_grp_dil)
    abs_log_dilution[i] = fabs(log_dilution[i]);
  std_loc = max(uID);
}
parameters{
  real<lower = 0> sigma_y;
  real<lower = 0> sigma_x;


  real<lower = 0> Bottom;
  real Span;
  real log_Inflec;
  real<lower = 0> Slope;

  real std_raw;

  vector[N_grp - 1] log_theta;
  vector[N_grp_dil] log_x;
}
model{
  vector[N_grp_dil] log_conc;

  sigma_y ~ normal(0, 10);
  sigma_x ~ normal(0, 1);
  Span ~ normal(16, 5);
  Bottom ~ normal(5, 1);
  log_Inflec ~ normal(inflec_mu, 1);
  Slope ~ normal(1, 0.3);
  std_raw ~ normal(0, 1);

  log_theta ~ normal(5, 10);


  for(i in 1:N_grp_dil){
    if(uID[i] == std_loc){
      log_conc[i] = (log(mu_Std + sigma_std * std_raw) + log_dilution[i]);
    } else {
      log_conc[i] = (log_theta[uID[i]] + log_dilution[i]);
    }
  }

  log_conc ~ normal(log_x, sigma_x);

  target += log(sigma_std) - log(sigma_std * std_raw + mu_Std);

  meas_Signal ~ normal(Bottom + Span * exp(-log1p_exp(-(log_x[dil_ID] - log_Inflec) * Slope)), sigma_y);

  zeroes ~ normal(Bottom, sigma_y);
}
