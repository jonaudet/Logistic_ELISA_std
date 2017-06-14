library(checkpoint)
checkpoint("2017-05-13", scanForPackages = F)
setSnapshot("2017-05-13")
if(!require(stringr)){
  install.packages("stringr")
  library(stringr)
}
if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}
if(!require(readxl)){
  install.packages("readxl")
  library(readxl)
}
if(!require(tidyr)){
  install.packages("tidyr")
  library(tidyr)
}
if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}
if(!require(rstan)){
  install.packages("rstan")
  library(rstan)
}
if(!require(rstudioapi)){
  install.packages("rstudioapi")
  library(rstudioapi)
}
if(!require(codetools)){
  install.packages("codetools")
  library(codetools)
}
if(!require(bayesplot)){
  install.packages("bayesplot")
  library(bayesplot)
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

unkn <- read_excel("E:/Cytokine/12Jun2017.xlsx", sheet = 2, skip = 1) %>%
  mutate(Sample = ifelse(str_detect(Sample, "S0"), "Std", Sample),
         uID = as.numeric(factor(Sample)),
         aID = as.numeric(factor(Assay))) %>%
  arrange(uID, Dilution) %>%
  mutate(Std = Sample == "Std",
         Dilution = 1 / Dilution)

unkn <- read_excel("E:/Cytokine/13Jun2017.xlsx", sheet = 2, skip = 1) %>%
  mutate(Sample = ifelse(str_detect(Sample, "S0"), "Std", Sample),
         uID = as.numeric(factor(Sample)),
         aID = as.numeric(factor(Assay))) %>%
  arrange(uID, Dilution) %>%
  mutate(Std = Sample == "Std",
         Dilution = 1 / Dilution) %>%
  bind_rows(unkn) %>%
  mutate(pID = as.numeric(factor(Plate_Name)))

unkn <- unkn %>%
  group_by(pID, aID) %>%
  mutate(Dilution = ifelse(is.na(Dilution), Concentration / max(Concentration, na.rm = T), Dilution)) %>%
  ungroup %>%
  group_by(Assay) %>%
  mutate(dID = as.numeric(factor(paste(Plate_Name, Sample, Dilution, sep = "_")))) %>%
  ungroup


ser_dilutions <- unkn %>%
  group_by(Assay) %>%
  mutate(dID = as.numeric(factor(paste(Plate_Name, Sample, Dilution, sep = "_")))) %>%
  distinct(dID, Dilution) %>%
  arrange(dID) %>%
  .$Dilution


# Plot of each plate's standard and corresponding unknowns
mutate(unkn, Dilution) %>%
  ggplot(aes(x = Concentration, y = Signal, colour = Plate_Name, group = uID)) +
  geom_point(alpha = 0.4) +
  #stat_summary(aes(fun.data = "mean"), geom = "line") +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  facet_wrap(~Assay, ncol = 3)


initial <- function(N_dil, N_plates, N_grp){
  inits <- list(std_raw = rnorm(1, 0, 1),
                mu_Bottom = abs(rnorm(1, 5, 1)),
                mu_Span = rnorm(1, 16, 5),
                mu_log_Inflec = rnorm(1, 5, 2),
                mu_log_Slope = abs(rnorm(1, 0, 0.15)),
                sigma_Bottom = abs(rnorm(1, 0, 1)),
                sigma_Span = abs(rnorm(1, 0, 1)),
                sigma_log_Inflec = abs(rnorm(1, 0, 1)),
                sigma_log_Slope = abs(rnorm(1, 0, 0.3)),
                log_theta = runif(N_grp - 1, -5, 6),
                log_x = runif(N_dil, -5, 6),
                sigma_x = rexp(1, 1),
                sigma_y = abs(rnorm(1, 0, 10)))
  return(inits)
}

# Run the model

mod <- stan_model("logistic_MSD_4p_JplateKassays.stan")

sep <- lapply(1:9, function(i){
  print(i)
  df <- filter(unkn, aID == i, Dilution != 0, Dilution != Inf) %>%
    mutate(uID = as.numeric(factor(Sample)),
           dID = as.numeric(factor(paste(Plate_Name, Sample, format(Dilution, scientific = F), sep = "_")))) %>%
    arrange(dID, uID, Dilution)
  inits <- lapply(1:4, function(x) initial(max(df$dID), 4, max(df$uID)))
  ser_dilutions <- df %>%
    mutate(dID = as.numeric(factor(paste(Plate_Name, Sample, format(Dilution, scientific = F), sep = "_")))) %>%
    distinct(dID, Dilution) %>%
    arrange(dID) %>%
    .$Dilution
  s_df <- df %>%
    distinct(dID, uID, .keep_all = T)
  mu_std <- max(df$Concentration, na.rm = T)
  zeroes <- filter(unkn, aID == i, Dilution == Inf) %>% .$Signal
  plate_zeroes <- filter(unkn, aID == i, Dilution == Inf) %>% .$pID
  inflec_mu <- mean(log(df$Concentration[df$Concentration > 0]), na.rm = T)
  res <- sampling(mod,
                  data = list(N = nrow(df),
                              J = max(df$pID),
                              N_grp = max(s_df$uID),
                              N_grp_dil = max(df$dID),
                              N_bot = length(zeroes),
                              dil_ID = df$dID,
                              uID = s_df$uID,
                              pID = s_df$pID,
                              meas_OD = log(df$Signal),
                              dilution = ser_dilutions,
                              mu_Std = mu_std,
                              sigma_std = 0.01 * mu_std,
                              zeroes = log(zeroes),
                              plate_zeroes = plate_zeroes,
                              inflec_mu = inflec_mu),
                  init = inits, chains = 4, sample_file = "dia",
                  iter = 3000, warmup = 500, refresh = 50, control = list(adapt_delta = 0.99, max_treedepth = 17))
  print(warnings())
  return(res)
})

## Check posterior
for(j in 1:9){
  print(j)
  log_x <- rstan::extract(sep[[j]], "log_x")$log_x

  ifn <- filter(unkn, aID == j, Dilution != 0, Dilution != Inf) %>%
    mutate(uID = as.numeric(factor(Sample)),
           dID = as.numeric(factor(paste(Plate_Name, Sample, format(Dilution, scientific = F), sep = "_")))) %>%
    arrange(dID, uID, Dilution)
  ifn$MedCalc <- 0
  ifn$Top <- 0
  ifn$Low <- 0

  for(i in 1:nrow(ifn)){
    if(ifn$dID[i] < 183){
      ifn$MedCalc[i] <- median(log_x[, ifn$dID[i]])
      ifn$Top[i] <- HDI(log_x[, ifn$dID[i]])$HigherHDI
      ifn$Low[i] <- HDI(log_x[, ifn$dID[i]])$LowerHDI
    }
  }


  print(filter(ifn, Std == T) %>%
    ggplot(aes(MedCalc, log(Signal), colour = Plate_Name)) +
    geom_point() +
    geom_segment(aes(yend = log(Signal), x = Low, xend = Top)) +
    labs(title = letters[j]))

}

log_theta <- array(0, dim = c(10000, 144, 9))

for(i in 1:9) log_theta[, , i] <- rstan::extract(sep[[i]], "log_theta")$log_theta

unkn$MedCalc <- 0
unkn$Low <- 0
unkn$Top <- 0

for(i in 1:nrow(unkn))
  if(unkn$uID[i] < 73){
    print(paste(unkn$uID[i], unkn$aID[i]))
    unkn$MedCalc[i] <- exp(median(log_theta[, unkn$uID[i], unkn$aID[i]]))
    unkn$Low[i] <- exp(HDI(log_theta[, unkn$uID[i], unkn$aID[i]])$LowerHDI)
    unkn$Top[i] <- exp(HDI(log_theta[, unkn$uID[i], unkn$aID[i]])$HigherHDI)
  }

unkn$Calc_Conc_Mean <- as.numeric(unkn$Calc_Conc_Mean)

filter(unkn, MedCalc != 0) %>%
  ggplot(aes(log(Calc_Conc_Mean), log(MedCalc), colour = Assay)) +
  geom_pointrange(aes(ymin = log(Low), ymax = log(Top)), alpha = 0.1) +
  geom_abline(slope = 1) +
  xlim(-9, 7) + ylim(-9, 7)

separate(unkn, Sample, c("Animal", "Day")) %>%
  mutate(Day = as.numeric(Day)) %>%
  filter(MedCalc != 0) %>%
  ggplot(aes(Day, log10(MedCalc), colour = Assay)) +
  geom_line() +
  geom_pointrange(aes(ymin = log10(Low), ymax = log10(Top))) +
  facet_wrap(~Animal) +
  coord_cartesian(ylim = c(-5, 5))

separate(unkn, Sample, c("Animal", "Day")) %>%
  mutate(Day = as.numeric(Day)) %>%
  filter(MedCalc != 0) %>%
  ggplot(aes(Day, log10(Calc_Conc_Mean), colour = Animal)) +
  geom_line() +
  geom_point() +
  facet_wrap(~Assay) +
  coord_cartesian(ylim = c(-5, 5))


## K Assays --------------------------------------------------------------------

mod <- stan_model("logistic_MSD_4p_JplateKassays.stan")

unkn <- read_excel("E:/Cytokine/12Jun2017.xlsx", sheet = 2, skip = 1) %>%
  mutate(Sample = ifelse(str_detect(Sample, "S0"), "Std", Sample),
         uID = as.numeric(factor(Sample)),
         aID = as.numeric(factor(Assay))) %>%
  arrange(uID, Dilution) %>%
  mutate(Std = Sample == "Std",
         Dilution = 1 / Dilution)

unkn <- read_excel("E:/Cytokine/13Jun2017.xlsx", sheet = 2, skip = 1) %>%
  mutate(Sample = ifelse(str_detect(Sample, "S0"), "Std", Sample),
         uID = as.numeric(factor(Sample)),
         aID = as.numeric(factor(Assay))) %>%
  arrange(uID, Dilution) %>%
  mutate(Std = Sample == "Std",
         Dilution = 1 / Dilution) %>%
  bind_rows(unkn) %>%
  mutate(pID = as.numeric(factor(Plate_Name)))

unkn <- unkn %>%
  group_by(pID, aID) %>%
  mutate(Dilution = ifelse(is.na(Dilution), Concentration / max(Concentration, na.rm = T), Dilution)) %>%
  ungroup %>%
  group_by(Assay) %>%
  mutate(dID = as.numeric(factor(paste(Plate_Name, Assay, Sample, Dilution, sep = "_")))) %>%
  ungroup


ser_dilutions <- unkn %>%
  group_by(Assay) %>%
  mutate(dID = as.numeric(factor(paste(Plate_Name, Assay, Sample, Dilution, sep = "_")))) %>%
  distinct(dID, Dilution) %>%
  arrange(dID) %>%
  .$Dilution

initial <- function(N_dil, N_plates, N_grp, N_assay){
  inits <- list(std_raw = rnorm(1, 0, 1),
                mu_Bottom = abs(rnorm(1, 5, 1)),
                mu_Span = rnorm(1, 16, 5),
                mu_log_Inflec = rnorm(1, 5, 2),
                mu_log_Slope = abs(rnorm(1, 0, 0.15)),
                sigma_plate_Bottom = abs(rnorm(N_assay, 0, 1)),
                sigma_plate_Span = abs(rnorm(N_assay, 0, 1)),
                sigma_plate_log_Inflec = abs(rnorm(N_assay, 0, 1)),
                sigma_plate_log_Slope = abs(rnorm(N_assay, 0, 0.3)),
                sigma_assay_Bottom = abs(rnorm(1, 0, 1)),
                sigma_assay_Span = abs(rnorm(1, 0, 1)),
                sigma_assay_log_Inflec = abs(rnorm(1, 0, 1)),
                sigma_assay_log_Slope = abs(rnorm(1, 0, 0.3)),
                log_theta = runif(N_grp - 1, -5, 6),
                log_x = runif(N_dil, -5, 6),
                sigma_x = rexp(1, 1),
                sigma_y = abs(rnorm(1, 0, 10)))
  return(inits)
}

df <- filter(unkn, Dilution != 0, Dilution != Inf)
inits <- lapply(1:4, function(x) initial(max(df$dID), 4, max(df$uID)))
s_df <- df %>%
  distinct(dID, uID, aID, .keep_all = T)
mu_std <- max(df$Concentration, na.rm = T)
zeroes <- filter(unkn, aID == i, Dilution == Inf) %>% .$Signal
plate_zeroes <- filter(unkn, aID == i, Dilution == Inf) %>% .$pID
inflec_mu <- mean(log(df$Concentration[df$Concentration > 0]), na.rm = T)
res <- sampling(mod,
                data = list(N = nrow(df),
                            J = max(df$pID),
                            N_grp = max(s_df$uID),
                            N_grp_dil = max(df$dID),
                            N_bot = length(zeroes),
                            dil_ID = df$dID,
                            uID = s_df$uID,
                            pID = s_df$pID,
                            meas_OD = log(df$Signal),
                            dilution = ser_dilutions,
                            mu_Std = mu_std,
                            sigma_std = 0.01 * mu_std,
                            zeroes = log(zeroes),
                            plate_zeroes = plate_zeroes,
                            inflec_mu = inflec_mu),
                init = inits, chains = 4, sample_file = "dia",
                iter = 3000, warmup = 500, refresh = 50, control = list(adapt_delta = 0.99, max_treedepth = 17))