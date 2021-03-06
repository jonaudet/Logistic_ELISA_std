data{
  int N;      //Total number of wells measured
  int N_grp;  //Total number of independent samples (1 sample has multiple dilutions which have multiple replicates)
  int N_grp_dil;  //Total number of dilutions and groups (so that N_grp_dil X # of replicates == N)
  int J;      //Total number of plates
  int N_bot; //Number of values at Conc = 0

  int<lower = 1, upper = N_grp> uID[N_grp_dil]; //Sample ID for each dilution
  int<lower = 1, upper = J> pID[N_grp_dil];     //Plate ID for each dilution
  int<lower = 1, upper = N_grp_dil> dil_ID[N];  //Dilution ID for each well
  vector<lower = 1>[N_grp_dil] Dil_order;

  vector[N] meas_OD;    //Optical density measured for each well
  vector[N_grp_dil] dilution;  //Dilution factor for each dilution

  real mu_Std;       //Concentration of the standard
  real<lower = 0> sigma_std;    //Uncertainty on the concentration of the standard

  vector[N_bot] zeroes;
  int plate_zeroes[N_bot];

  real inflec_mu;
  
  int N_mix;
  vector[N_mix] means;
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
  //real<lower = 0> sigma_x;


  real mu_Bottom;
  real mu_Span;
  real mu_log_Inflec;
  real mu_log_Slope;

  real<lower = 0> sigma_Bottom;
  real<lower = 0> sigma_Span;
  real<lower = 0> sigma_log_Inflec;
  real<lower = 0> sigma_log_Slope;

  vector[J] raw_Bottom;
  vector[J] raw_Span;
  vector[J] raw_log_Inflec;
  vector[J] raw_log_Slope;

  real std_raw;

  vector[N_grp - 1] log_theta;
  vector[N_grp_dil] log_x;
}
model{
  vector[N_grp_dil] log_conc;
  vector[J] Bottom;
  vector[J] Span;
  vector[J] log_Inflec;
  vector[J] Slope;

  for(i in 1:N_grp_dil){
    if(uID[i] == std_loc){
      log_conc[i] = (log(mu_Std + sigma_std * std_raw) + log_dilution[i]);
    } else {
      log_conc[i] = (log_theta[uID[i]] + log_dilution[i]);
    }
  }

  Bottom = mu_Bottom + sigma_Bottom * raw_Bottom;
  Span = mu_Span + sigma_Span * raw_Span;
  log_Inflec = mu_log_Inflec + sigma_log_Inflec * raw_log_Inflec;
  Slope = exp(mu_log_Slope + sigma_log_Slope * raw_log_Slope);

  sigma_y ~ normal(0, 1);
  //sigma_x ~ normal(0, 1);
  mu_Span ~ normal(10, 2);
  mu_Bottom ~ normal(5, 1);
  mu_log_Inflec ~ normal(inflec_mu, 2);
  mu_log_Slope ~ normal(0, 0.12);

  sigma_Bottom ~ normal(0, 1);
  sigma_Span ~ normal(0, 1);
  sigma_log_Inflec ~ normal(0, 1);
  sigma_log_Slope ~ normal(0, 0.3);

  raw_Bottom ~ normal(0, 1);
  raw_Span ~ normal(0, 1);
  raw_log_Inflec ~ normal(0, 1);
  raw_log_Slope ~ normal(0, 1);

  std_raw ~ normal(0, 1);

  //log_theta prior, flat-ish over -4 to 11
  
  for(j in 1:(N_grp - 1)){
    vector[N_mix] probs;
    for(i in 1:N_mix){
      probs[i] = log(1.0 / N_mix) + normal_lpdf(log_theta[j] | means[i], 1);
    }
    
    target += log_sum_exp(probs);
  }

  log_conc ~ normal(log_x, 0.05 + 0.01 * Dil_order);

  target += log(sigma_std) - log(sigma_std * std_raw + mu_Std);

  meas_OD ~ normal(Bottom[pID[dil_ID]] + Span[pID[dil_ID]] .* exp(-log1p_exp(-(log_x[dil_ID] - log_Inflec[pID[dil_ID]]) .* Slope[pID[dil_ID]])), sigma_y);

  zeroes ~ normal(Bottom[plate_zeroes], sigma_y);
}
