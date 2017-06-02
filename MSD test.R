if(!require(checkpoint)){
  install.packages("checkpoint")
  library(checkpoint)
}
checkpoint("2017-05-13", scanForPackages = F)
if(!(require(readr) & require(ggplot2) & require(tidyr) & require(plyr) & require(dplyr) & require(rstan) & require(rstudioapi) & require(codetools) & require(readxl))){
  install.packages(c("rstan", "ggplot2", "tidyr", "dplyr", "plyr", "readr", "rstudioapi", "codetools", "readxl"))

  library(rstan)
  library(readr)
  library(ggplot2)
  library(plyr)
  library(tidyr)
  library(dplyr)
  library(rstudioapi)
  library(codetools)
  library(readxl)
}
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

HDI <- function(Values, Interval = 0.95){
  Values <- Values[!is.na(Values)]
  intSize <- as.integer(length(Values) * Interval)
  startMax <- as.integer(length(Values) - (length(Values) * Interval))
  ordered <- Values[sort.list(Values)]
  low <- 1
  diff <- Inf
  for(i in 1:startMax)
    if(ordered[i + intSize] - ordered[i] < diff){
      low <- i
      diff <- ordered[i + intSize] - ordered[i]
    }
  return(data.frame(LowerHDI = ordered[low], HigherHDI = ordered[low + intSize]))
}

unkn <- read_excel("G:/Cytokine16May.xlsx", sheet = 2) %>%
  mutate(Sample = ifelse(Sample %in% paste("Std", 1:8, sep = ""), "Std", Sample),
         uID = as.numeric(factor(Sample)),
         aID = as.numeric(factor(Assay))) %>%
  arrange(uID, Dilution) %>%
  mutate(Std = Sample == "Std",
         Initial = Dilution * Concentration,
         Dilution = 1 / Dilution) %>%
  group_by(Assay) %>%
  mutate(dID = as.numeric(factor(paste(Sample, Dilution, sep = "_")))) %>%
  ungroup

ser_dilutions <- unkn %>%
  group_by(Assay) %>%
  mutate(dID = as.numeric(factor(paste(Sample, Dilution, sep = "_")))) %>%
  distinct(dID, Dilution) %>%
  arrange(dID) %>%
  .$Dilution

unkn <- group_by(unkn, Assay) %>%
  mutate(dID = as.numeric(factor(paste(Sample, Dilution, sep = "_")))) %>%
  ungroup

# Plot of each plate's standard and corresponding unknowns
mutate(unkn, Dilution) %>%
  ggplot(aes(x = Concentration, y = Signal, colour = Std, group = uID)) +
  geom_point(alpha = 0.4) +
  #stat_summary(aes(fun.data = "mean"), geom = "line") +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  facet_wrap(~Assay, ncol = 3)

mutate(unkn, Conc = 4500 * Dilution) %>%
  filter(Std == TRUE) %>%
  ggplot(aes(x = Conc, y = Signal, colour = Assay)) +
  geom_point(alpha = 0.4) +
  stat_summary(aes(fun.data = "mean"), geom = "line") +
  scale_x_log10() +
  theme_bw()


initial <- function(N_dil, N_plates, N_grp){
  inits <- list(std_raw = rnorm(1, 0, 1),
                sigma_y = abs(rnorm(1, 0, 10)),
                Bottom = abs(rnorm(N_plates, 5, 1)),
                Span = rnorm(N_plates, 16, 5),
                log_Inflec = rnorm(N_plates, 5, 2),
                Slope = abs(rnorm(N_plates, 1, 0.3)),
                log_theta = runif(N_grp - 1, -5, 6),
                sigma_x = rexp(1, 1),
                sigma_y = abs(rnorm(1, 0, 10)))
  return(inits)
}

# Run the model

mod <- stan_model("logistic_MSD_4p_1plate.stan")

sep <- lapply(1:9, function(i){
  print(i)
  df <- filter(unkn, aID == i, !is.na(Dilution)) %>%
    mutate(uID = as.numeric(factor(Sample)),
           dID = as.numeric(factor(paste(Sample, format(Dilution, scientific = F), sep = "_")))) %>%
    arrange(dID, uID, Dilution)
  inits <- lapply(1:4, function(x) initial(max(df$dID), 1, max(df$uID)))
  ser_dilutions <- df %>%
    mutate(dID = as.numeric(factor(paste(Sample, format(Dilution, scientific = F), sep = "_")))) %>%
    distinct(dID, Dilution) %>%
    arrange(dID) %>%
    .$Dilution
  s_df <- df %>%
    distinct(dID, uID)
  mu_std <- median(df$Initial, na.rm = T)
  zeroes <- filter(unkn, aID == i, is.na(Dilution)) %>% .$Signal
  inflec_mu <- mean(df$Concentration, na.rm = T)
  res <- sampling(mod,
                  data = list(N = nrow(df),
                              N_grp = max(s_df$uID),
                              N_grp_dil = max(df$dID),
                              N_bot = length(zeroes),
                              dil_ID = df$dID,
                              uID = s_df$uID,
                              meas_Signal = log(df$Signal),
                              dilution = ser_dilutions,
                              mu_Std = mu_std,
                              sigma_std = 0.01 * mu_std,
                              zeroes = log(zeroes),
                              inflec_mu = log(inflec_mu)),
                  init = inits, chains = 4, sample_file = "dia",
                  iter = 3000, warmup = 500, refresh = 50, control = list(adapt_delta = 0.99, max_treedepth = 17))
  print(warnings())
  return(res)
})

logx <- rstan::extract(res, "log_x")$log_x
df$logx <- 0
df$logx_low <- 0
df$logx_top <- 0
for(i in 1:nrow(df)){
  df$logx[i] <- median(logx[, df$dID[i]])
  df$logx_low[i] <- HDI(logx[, df$dID[i]], Interval = 0.97)$LowerHDI
  df$logx_top[i] <- HDI(logx[, df$dID[i]], Interval = 0.97)$HigherHDI
}
ggplot(df, aes(logx, Signal)) +
  geom_segment(aes(yend = Signal, x = logx_low, xend = logx_top, colour = Std)) +
  geom_point(aes(colour = Std)) +
  scale_y_log10() +
  geom_hline(yintercept = 250.5449, linetype = 3)

log_theta <- array(0, dim = c(10000, 36, 9))

for(i in 1:9) log_theta[, , i] <- rstan::extract(sep[[i]], "log_theta")$log_theta

unkn$MedCalc <- 0

for(i in 1:nrow(unkn))
  if(unkn$uID[i] < 37)
    #print(paste(unkn$uID[i], unkn$aID[i]))
    unkn$MedCalc[i] <- exp(median(log_theta[, unkn$uID[i], unkn$aID[i]]))

unkn$Calc_Conc_Mean <- as.numeric(unkn$Calc_Conc_Mean)

filter(unkn, MedCalc != 0) %>%
ggplot(aes(log(Calc_Conc_Mean), log(MedCalc))) + geom_point(alpha = 0.1) + geom_abline(slope = 1) + xlim(-9, 7) + ylim(-9, 7)

