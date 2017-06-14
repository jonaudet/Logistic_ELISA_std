data{
  int N;      //Total number of wells measured
  int N_grp;  //Total number of independent samples (1 sample has multiple dilutions which have multiple replicates)
  int N_grp_dil;  //Total number of dilutions and groups (so that N_grp_dil X # of replicates == N)
  int J;      //Total number of plates
  int N_bot; //Number of values at Conc = 0
  int K;     //Total number of assays

  int<lower = 1, upper = N_grp> uID[N_grp_dil]; //Sample ID for each dilution
  int<lower = 1, upper = J> pID[N_grp_dil];     //Plate ID for each dilution
  int<lower = 1, upper = N_grp_dil> dil_ID[N];  //Dilution ID for each well
  int<lower = 1, upper = K> assay_ID[N_grp_dil];//Assay ID for each dilution

  vector[N] meas_OD;    //Optical density measured for each well
  vector[N_grp_dil] dilution;  //Dilution factor for each dilution

  real mu_Std;       //Concentration of the standard
  real<lower = 0> sigma_std;    //Uncertainty on the concentration of the standard

  vector[N_bot] zeroes;
  int plate_zeroes[N_bot];
  int assay_zeroes[N_bot];

  real inflec_mu;
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
  real<lower = 0> sigma_y;  //Variation in measurement
  real<lower = 0> sigma_x;  //Variation in dilution


  real mu_Bottom;  //Overall average Bottom
  real mu_Span;    //Overall average Span
  real mu_log_Inflec;  //Overall average inflection point
  real mu_log_Slope;  //Overall average slope

  row_vector<lower = 0>[K] sigma_plate_Bottom;       //Between-plate variation in Bottom for each assay
  row_vector<lower = 0>[K] sigma_plate_Span;       //Between-plate variation in Span for each assay
  row_vector<lower = 0>[K] sigma_plate_log_Inflec;       //Between-plate variation in Inflection point for each assay
  row_vector<lower = 0>[K] sigma_plate_log_Slope;       //Between-plate variation in Slope for each assay

  //Z-value for each plate within each assay
  matrix[J, K] raw_plate_Bottom;
  matrix[J, K] raw_plate_Span;
  matrix[J, K] raw_plate_log_Inflec;
  matrix[J, K] raw_plate_log_Slope;

  //Between-assay variation
  real<lower = 0> sigma_assay_Bottom;
  real<lower = 0> sigma_assay_Span;
  real<lower = 0> sigma_assay_log_Inflec;
  real<lower = 0> sigma_assay_log_Slope;

  //Z-value for each assay
  row_vector[K] raw_assay_Bottom;
  row_vector[K] raw_assay_Span;
  row_vector[K] raw_assay_log_Inflec;
  row_vector[K] raw_assay_log_Slope;

  real std_raw;

  vector[N_grp - 1] log_theta;
  vector[N_grp_dil] log_x;
}
model{
  vector[N_grp_dil] log_conc;
  matrix[J, K] Bottom;
  matrix[J, K] Span;
  matrix[J, K] log_Inflec;
  matrix[J, K] Slope;
  vector[N] calc_ODs;

  for(i in 1:N_grp_dil){
    if(uID[i] == std_loc){
      log_conc[i] = (log(mu_Std + sigma_std * std_raw) + log_dilution[i]);
    } else {
      log_conc[i] = (log_theta[uID[i]] + log_dilution[i]);
    }
  }

  Bottom = rep_matrix(mu_Bottom, J, K) +
            rep_matrix(sigma_plate_Bottom, J) .* raw_plate_Bottom +
            rep_matrix(sigma_assay_Bottom, J, K) .* rep_matrix(raw_assay_Bottom, J);
  Span = rep_matrix(mu_Span, J, K) +
            rep_matrix(sigma_plate_Span, J) .* raw_plate_Span +
            rep_matrix(sigma_assay_Span, J, K) .* rep_matrix(raw_assay_Span, J);
  log_Inflec = rep_matrix(mu_log_Inflec, J, K) +
            rep_matrix(sigma_plate_log_Inflec, J) .* raw_plate_log_Inflec +
            rep_matrix(sigma_assay_log_Inflec, J, K) .* rep_matrix(raw_assay_log_Inflec, J);
  Slope = exp(rep_matrix(mu_log_Slope, J, K) +
            rep_matrix(sigma_plate_log_Slope, J) .* raw_plate_log_Slope +
            rep_matrix(sigma_assay_log_Slope, J, K) .* rep_matrix(raw_assay_log_Slope, J));

  sigma_y ~ normal(0, 1);
  sigma_x ~ normal(0, 1);
  mu_Span ~ normal(10, 2);
  mu_Bottom ~ normal(5, 1);
  mu_log_Inflec ~ normal(inflec_mu, 2);
  mu_log_Slope ~ normal(0, 0.12);

  sigma_plate_Bottom ~ normal(0, 1);
  sigma_plate_Span ~ normal(0, 1);
  sigma_plate_log_Inflec ~ normal(0, 1);
  sigma_plate_log_Slope ~ normal(0, 0.3);

  to_vector(raw_plate_Bottom) ~ normal(0, 1);
  to_vector(raw_plate_Span) ~ normal(0, 1);
  to_vector(raw_plate_log_Inflec) ~ normal(0, 1);
  to_vector(raw_plate_log_Slope) ~ normal(0, 1);

  sigma_assay_Bottom ~ normal(0, 1);
  sigma_assay_Span ~ normal(0, 1);
  sigma_assay_log_Inflec ~ normal(0, 1);
  sigma_assay_log_Slope ~ normal(0, 0.3);

  raw_assay_Bottom ~ normal(0, 1);
  raw_assay_Span ~ normal(0, 1);
  raw_assay_log_Inflec ~ normal(0, 1);
  raw_assay_log_Slope ~ normal(0, 1);

  std_raw ~ normal(0, 1);

  log_theta ~ normal(5, 10);

  log_conc ~ normal(log_x, sigma_x);

  target += log(sigma_std) - log(sigma_std * std_raw + mu_Std);

  for(i in 1:N)
    calc_ODs[i] = Bottom[pID[dil_ID[i]], assay_ID[dil_ID[i]]] +
                                          Span[pID[dil_ID[i]], assay_ID[dil_ID[i]]] .* exp(-log1p_exp(-(log_x[dil_ID[i]] - log_Inflec[pID[dil_ID[i]], assay_ID[dil_ID[i]]]) .* Slope[pID[dil_ID[i]], assay_ID[dil_ID[i]]]));

  meas_OD ~ normal(calc_ODs, sigma_y);

  zeroes ~ normal(to_vector(Bottom[plate_zeroes, assay_zeroes]), sigma_y);
}
