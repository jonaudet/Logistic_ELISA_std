data{
  int N_unkn;           //Total number of unknown values
  int N_unkn_grp;       //Total number of unknown samples (should be N_unkn / 4)
  int N_plates;
  int uID[N_unkn];      //Sample number for each unknown value
  vector[N_unkn] measured_OD; //Optical Density of each unknown or standard
  vector[N_unkn] ser_dilutions; //The serial dilutions from the start dilution for each unknown;
  real mu_Std;
  real<lower = 0> sigma_std;
  int pID[N_unkn]; //plate # for each unknown (to match it to the proper standard)
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
  real pred_std_raw;
  real<lower = 0> sigma_unkn;
  
    real mu_log_Bottom;
  real<lower = 0> sigma_log_Bottom;
  
  real<lower = mu_log_Bottom> mu_log_Span;
  real<lower = 0> sigma_log_Span;
  
  // real mu_log_Inflec;
  real<lower = 0> sigma_log_Inflec;
  
  real mu_log_Slope;
  real<lower = 0> sigma_log_Slope;
  // 
  // real mu_log_Asym;
  real<lower = 0> sigma_log_Asym;

  vector[N_unkn_grp - 1] log_theta;  // log of each unknown's predicted concentration
  //vector<lower = 0>[N_unkn_grp] sigma_x;
  real<lower = 0> sigma_x;
  
  vector[N_plates] z_log_Asyms;
  vector[N_plates] z_log_Spans;
  vector[N_plates] z_log_Bottoms;
  vector[N_plates] z_log_Inflecs;
  vector[N_plates] z_log_Slopes;
  vector[N_unkn] log_x_raw; // log of each unknown's concentration
  cholesky_factor_corr[2] L;
  vector<lower =  0>[2] L_sigma;
  vector[2] alpha;
  vector[2] mu;
}
transformed parameters{

  vector<lower = 0>[N_plates] Asyms;
  vector<lower = 0>[N_plates] Spans;
  vector<lower = 0>[N_plates] Bottoms;
  vector[N_plates] log_Inflecs;
  vector<lower = 0>[N_plates] Slopes;
  real<upper = 12> mu_log_Asym;
  real mu_log_Inflec;
  
  {
    vector[2] temp;
    
    temp <- mu + diag_pre_multiply(L_sigma, L) * alpha;
    
    mu_log_Asym <- temp[1];
    mu_log_Inflec <- temp[2];
  }
  
  
  Bottoms <- exp(mu_log_Bottom + z_log_Bottoms * sigma_log_Bottom);
  Spans <- exp(mu_log_Span + z_log_Spans * sigma_log_Span);
  log_Inflecs <- mu_log_Inflec + z_log_Inflecs * sigma_log_Inflec;
  Asyms <- exp(mu_log_Asym + z_log_Asyms * sigma_log_Asym);
  Slopes <- exp(mu_log_Slope + z_log_Slopes * sigma_log_Slope);
  

  #print(pred_std);
  #print(log_Undil);
      
}
model{
  vector[N_unkn] unkn_cOD;
  vector[N_unkn] log_x;
  vector[N_unkn] log_Undil;
  real pred_std;
  
  
  pred_std <- mu_Std + pred_std_raw * sigma_std;
  
  for(i in 1:N_unkn){
    if(uID[i] == std_loc){
      log_Undil[i] <- log(pred_std);
    } else {
      log_Undil[i] <- log_theta[uID[i]];
    }
  }
  
  log_x <- (log_ser_dilutions + log_Undil) + (sigma_x * abs_log_ser_dilutions) .* log_x_raw;

  for(i in 1:N_unkn){
    // unkn_cOD[i] <- Bottoms[pID[i]] + (Spans[pID[i]] / (1 + exp((log_Inflecs[pID[i]] - log_x[i]) * Slopes[pID[i]])) ^ Asyms[pID[i]]);
    unkn_cOD[i] <- Bottoms[pID[i]] + (Spans[pID[i]] / (1 + exp((log_Inflecs[pID[i]] - log_x[i]) * Slopes[pID[i]])));
  }
  
  z_log_Bottoms ~ normal(0, 1);
  mu_log_Bottom ~ normal(log(0.05), 2);
  sigma_log_Bottom ~ normal(0, 1);
  
  z_log_Spans ~ normal(0, 1);
  mu_log_Span ~ normal(1.2, 0.1);
  sigma_log_Span ~ normal(0, 1);
  
  z_log_Inflecs ~ normal(0, 1);
  // mu_log_Inflec ~ normal(0, 3);
  mu[2] ~ normal(0, 1);
  sigma_log_Inflec ~ normal(0, 1);
  
  z_log_Slopes ~ normal(0, 1);
  mu_log_Slope ~ normal(0, 1);
  sigma_log_Slope ~ normal(0, 1);
  
  z_log_Asyms ~ normal(0, 1);
  // mu_log_Asym ~ normal(0, 2);
  mu[1] ~ normal(0, 1);
  sigma_log_Asym ~ normal(0, 1);
  
  log_theta ~ uniform(-10, 15);
  sigma_x ~ normal(0, 1);

  alpha ~ normal(0, 1);
  L ~ lkj_corr_cholesky(4);
  L_sigma ~ normal(0, 1);

  sigma_unkn ~ normal(0, 1);
  
  
  //Multilevel unknown estimation
  //log_x ~ normal(log_ser_dilutions + log_Undil, sigma_x * abs_log_ser_dilutions);
  log_x_raw ~ normal(0, 1);
  
  pred_std_raw ~ normal(0, 1);

  measured_OD ~ normal(unkn_cOD, sigma_unkn);
}
generated quantities{
  vector[N_unkn_grp - 1] theta;
  vector[N_unkn] x;
  
  theta <- exp(log_theta);
  x <- exp(log_x);
}