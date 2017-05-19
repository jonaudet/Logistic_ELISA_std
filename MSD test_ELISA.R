if(!require(checkpoint)){
  install.packages("checkpoint")
  library(checkpoint)
}
checkpoint("2017-05-13", scanForPackages = F)
setSnapshot("2017-05-13")
if(!(require(readr) & require(ggplot2) & require(tidyr) & require(plyr) & require(dplyr) & require(rstan) & require(rstudioapi) & require(codetools) & require(readxl) & require(bayesplot))){
  install.packages(c("rstan", "ggplot2", "tidyr", "dplyr", "plyr", "readr", "rstudioapi", "codetools", "readxl", "bayesplot"))

  library(rstan)
  library(readr)
  library(ggplot2)
  library(plyr)
  library(tidyr)
  library(dplyr)
  library(rstudioapi)
  library(codetools)
  library(readxl)
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

LoadData <- function(Folder = "../Data/qELISA-2", DataFile = "Plates.xlsx", DilFile = "Dilutions.xlsx", LayoutFile = "Layout.xlsx"){
  DataPath <- paste(Folder, DataFile, sep = "/")
  DilPath <- paste(Folder, DilFile, sep = "/")
  LayoutPath <- paste(Folder, LayoutFile, sep = "/")

  Data <- data_frame()
  for(i in excel_sheets(DataPath))
    Data <- read_excel(DataPath, sheet = i) %>%
    gather("Column", "OD", 2:length(.)) %>%
    arrange(Row, Column) %>%
    dplyr::mutate(Plate = i, Rep = rep(c("OD1", "OD2"), length(.$Column)/2)) %>%
    unite(Well, Plate, Row, Column, sep = "_") %>%
    bind_rows(Data, .)

  Dilutions <- data_frame()
  for(i in excel_sheets(DilPath))
    Dilutions <- read_excel(DilPath, sheet = i) %>%
    gather("Column", "Dilution", 2:length(.)) %>%
    arrange(Row, Column) %>%
    dplyr::mutate(Plate = i, Rep = rep(c("OD1", "OD2"), length(.$Column)/2)) %>%
    unite(Well, Plate, Row, Column, sep = "_") %>%
    filter(!is.na(Dilution)) %>%
    mutate(Dilution = 1 / Dilution) %>%
    bind_rows(Dilutions, .)

  Samp <- data_frame()
  for(i in excel_sheets(LayoutPath))
    Samp <- read_excel(LayoutPath, sheet = i) %>%
    gather("Column", "Samp", 2:length(.)) %>%
    arrange(Row, Column) %>%
    dplyr::mutate(Plate = i, Rep = rep(c("OD1", "OD2"), length(.$Column)/2)) %>%
    unite(Well, Plate, Row, Column, sep = "_") %>%
    filter(!is.na(Samp)) %>%
    bind_rows(Samp, .)

  Data <- Data %>%
    left_join(Dilutions) %>%
    left_join(Samp) %>%
    separate(Well, c("Plate", "Well"), sep = "_", extra = "merge") %>%
    filter(!is.na(Samp)) %>%
    dplyr::mutate(Well = gsub("_", "", Well))

  return(Data)
}

unkn <- LoadData(Folder = "D:", DataFile = "Signal.xlsx", DilFile = "Dilution_IgG.xlsx", LayoutFile = "Layout_IgG.xlsx") %>%
  mutate(uID = as.numeric(factor(Samp)),
         pID = as.numeric(factor(Plate)),
         dID = as.numeric(factor(paste(Plate, Samp, format(Dilution, scientific = F), sep = "_")))) %>%
  arrange(dID, Dilution) %>%
  mutate(Std = Samp == "std")

ser_dilutions <- unkn %>%
  mutate(dID = as.numeric(factor(paste(Plate, Samp, Dilution, sep = "_")))) %>%
  distinct(dID, Dilution) %>%
  arrange(dID) %>%
  .$Dilution

unkn <- mutate(unkn, dID = as.numeric(factor(paste(Plate, Samp, Dilution, sep = "_"))))

mutate(unkn, Conc = 2040 * Dilution) %>%
  ggplot(aes(x = Conc, y = log(OD), colour = Std, group = uID)) +
  geom_point(alpha = 0.4) +
  stat_summary(aes(fun.data = "mean"), geom = "line") +
  scale_x_log10() +
  theme_bw() +
  facet_wrap(~Plate, ncol = 3)

mutate(unkn, Conc = 4500 * Dilution) %>%
  filter(Std == TRUE) %>%
  ggplot(aes(x = Conc, y = OD, colour = Plate)) +
  geom_point(alpha = 0.4) +
  stat_summary(aes(fun.data = "mean"), geom = "line") +
  scale_x_log10() +
  theme_bw()


initial <- function(N_dil, N_plates, N_grp){
  inits <- list(std_raw = rnorm(1, 0, 1),
                sigma_y = abs(rnorm(1, 0, 1)),
                Bottom = abs(rnorm(N_plates, 5, 1)),
                Span = rnorm(N_plates, 8, 1),
                log_Inflec = rnorm(N_plates, 0, 1),
                Slope = abs(rnorm(N_plates, 1, 0.5)),
                log_theta = runif(N_grp - 1, -5, 6),
                sigma_x = rexp(1, 1),
                sigma_OD = abs(rnorm(1, 0, 0.2)))
  return(inits)
}

# Run the model

mod <- stan_model("logistic_MSD_ELISA_4p_1plate.stan")

sep <- lapply(1:2, function(i){
  print(i)
  df <- filter(unkn, pID == i) %>%
    mutate(uID = as.numeric(factor(Samp)),
           pID = as.numeric(factor(Plate)),
           dID = as.numeric(factor(paste(Samp, format(Dilution, scientific = F), sep = "_")))) %>%
    arrange(dID, uID, Dilution) %>%
    mutate(Std = Samp == "std")
  inits <- lapply(1:4, function(x) initial(max(df$dID), 1, max(df$uID)))
  ser_dilutions <- df %>%
    mutate(dID = as.numeric(factor(paste(Samp, format(Dilution, scientific = F), sep = "_")))) %>%
    distinct(dID, Dilution) %>%
    arrange(dID) %>%
    .$Dilution
  s_df <- df %>%
    distinct(dID, uID)
  res <- sampling(mod,
                  data = list(N = nrow(df),
                              N_grp = max(s_df$uID),
                              N_grp_dil = max(df$dID),
                              dil_ID = df$dID,
                              uID = s_df$uID,
                              meas_OD = log(df$OD),
                              dilution = ser_dilutions,
                              mu_Std = 4500,
                              sigma_std = 200),
                  init = inits, chains = 4, sample_file = "dia",
                  iter = 2000, warmup = 500, refresh = 50, control = list(adapt_delta = 0.99, max_treedepth = 14))
  return(res)
})


check_divergent <- function(stan_res){
  np <- nuts_params(stan_res)
  divergences <- filter(np, Parameter == "divergent__", Iteration > 500)
  number <- sum(divergences$Value)
  return(number)
}

Divergences <- sapply(1:2, function(x) check_divergent(sep[[x]]))


output <- unkn
out_sep <- bind_rows(lapply(1:2, function(i){
  cat(i)
  df <- filter(unkn, pID == i) %>%
    mutate(uID = as.numeric(factor(Samp)),
           pID = as.numeric(factor(Plate)),
           dID = as.numeric(factor(paste(Samp, format(Dilution, scientific = F), sep = "_")))) %>%
    arrange(dID, uID, Dilution) %>%
    mutate(Std = Samp == "std")

  cOD <- exp(rstan::extract(sep[[i]], "log_x")$log_x)
  df$Median <- rep(apply(cOD, 2, median), each = 2)
  errors <- bind_rows(apply(cOD, 2, HDI))
  df$TopHDI <- rep(errors$HigherHDI, each = 2)
  df$LowHDI <- rep(errors$LowerHDI, each = 2)
  df$Conc <- df$Median
  return(df)
}))

out_sep <- mutate(out_sep, uID = as.numeric(factor(Samp)))
ggplot(out_sep, aes(Conc, log(OD))) +
  scale_x_log10(breaks = 10^seq(floor(log10(min(out_sep$Conc))), ceiling(log10(max(out_sep$Conc))), by = 1)) +
  coord_cartesian(xlim = c(1e-4, 400)) +
  #geom_point(aes(colour = factor(uID))) +
  geom_text(aes(label = Samp, colour = factor(uID))) +
  geom_errorbarh(aes(xmin = LowHDI, xmax = TopHDI, colour = factor(uID))) +
  scale_colour_discrete(guide = "none") +
  facet_wrap(~Plate, ncol = 3)

out_sep <- bind_rows(lapply(1:2, function(i){
  print(i)
  out <- unkn %>%
    filter(Samp != "std", pID == i) %>%
    mutate(uID = as.numeric(factor(Samp)),
           pID = as.numeric(factor(Plate))) %>%
    arrange(uID) %>%
    distinct(uID, .keep_all = T) %>%
    #separate(Samp, c("Group", "Samp"), sep = "-") %>%
    separate(Samp, c("Unit", "Week"), sep = "_", remove = F) %>%
    mutate(Week = as.numeric(Week))

  theta <- exp(rstan::extract(sep[[i]], "log_theta")$log_theta)
  out$Conc <- apply(theta, 2, median)[out$uID]
  errors <- bind_rows(apply(theta, 2, HDI))
  out$TopHDI <- errors$HigherHDI[out$uID]
  out$LowHDI <- errors$LowerHDI[out$uID]
  return(out)
}))

out_sep <- mutate(out_sep, pID = as.numeric(factor(Plate)),
                  Method = "Individual",
                  uID = as.numeric(factor(Samp)))
ggplot(out_sep, aes(Week, Conc, colour = Unit, fill = Unit)) +
  geom_pointrange(aes(ymin = LowHDI, ymax = TopHDI, shape = Unit)) +
  geom_line() +
  scale_y_log10(breaks = 10^seq(-12, 6)) +
  annotation_logticks(sides = "l") +
  coord_cartesian(ylim = c(0.001, 5e4)) +
  xlim(0, NA) +
  facet_wrap(~Plate) +
  theme_bw()

# Multiplate -------------------------------------------------------------------

ser_dilutions <- unkn %>%
  mutate(dID = as.numeric(factor(paste(Plate, Samp, format(Dilution, scientific = F), sep = "_")))) %>%
  distinct(dID, Dilution) %>%
  arrange(dID) %>%
  .$Dilution

dil_unkn <- unkn %>%
  distinct(uID, pID, dID, Dilution)

mod <- stan_model("logistic_X_4p_Jplate.stan")

res2 <- sampling(mod,
                 data = list(N = nrow(unkn),
                             N_grp = max(unkn$uID),
                             N_grp_dil = max(unkn$dID),
                             uID = dil_unkn$uID,
                             dil_ID = unkn$dID,
                             meas_OD = unkn$OD,
                             dilution = ser_dilutions,
                             mu_Std = 2000,
                             sigma_std = 5,
                             J = max(unkn$pID),
                             pID = dil_unkn$pID),
                 init = inits, chains = 4,#diagnostic_file = "dia",
                 iter = 2000, warmup = 500, refresh = 50, control = list(adapt_delta = 0.95))

stan_ess(res2)
stan_ac(res2)
stan_rhat(res2)

# Look at the curves and the unknowns, the standards are shown in red all other colors are unknowns

output <- unkn

cOD <- exp(rstan::extract(res2, "log_x")$log_x)
output$Median <- apply(cOD, 2, median)[output$dID]
errors <- bind_rows(apply(cOD, 2, HDI))
output$TopHDI <- errors$HigherHDI[output$dID]
output$LowHDI <- errors$LowerHDI[output$dID]
output$Conc <- output$Median

ggplot(filter(output, Std == T), aes(Conc, OD)) +
  scale_x_log10(breaks = 10^seq(floor(log10(min(output$Conc))), ceiling(log10(max(output$Conc))), by = 1)) +
  coord_cartesian(xlim = c(1e-4, 35), ylim = c(0, 4)) +
  #geom_point(aes(colour = factor(uID))) +
  geom_text(aes(label = Samp, colour = factor(uID))) +
  geom_errorbarh(aes(xmin = LowHDI, xmax = TopHDI, colour = factor(uID))) +
  scale_colour_discrete(guide = "none") +
  facet_wrap(~Plate, ncol = 3)

ggplot(output, aes(Conc, OD)) +
  scale_x_log10(breaks = 10^seq(floor(log10(min(output$Conc))), ceiling(log10(max(output$Conc))), by = 1)) +
  coord_cartesian(xlim = c(1e-4, 35), ylim = c(0, 4)) +
  #geom_point(aes(colour = factor(uID))) +
  geom_text(aes(label = Samp, colour = factor(uID))) +
  geom_errorbarh(aes(xmin = LowHDI, xmax = TopHDI, colour = factor(uID))) +
  scale_colour_discrete(guide = "none") +
  facet_wrap(~Plate, ncol = 3)

# Plotting the output data (theta) as it is meant to

out <- unkn %>%
  filter(Samp != "std") %>%
  arrange(uID) %>%
  group_by(uID) %>%
  top_n(1, OD) %>%
  ungroup %>%
  #separate(Samp, c("Group", "Samp"), sep = "-") %>%
  separate(Samp, c("Unit", "Day"), sep = "_") %>%
  mutate(Day = as.numeric(Day))

theta <- exp(rstan::extract(res2, "log_theta")$log_theta)
out$Conc <- apply(theta, 2, median)[out$uID]
errors <- bind_rows(apply(theta, 2, HDI))
out$TopHDI <- errors$HigherHDI[out$uID]
out$LowHDI <- errors$LowerHDI[out$uID]

out$Method <- "Mulitlevel"

ggplot(out, aes(Day, Conc, colour = Unit, fill = Unit)) +
  geom_pointrange(aes(ymin = LowHDI, ymax = TopHDI, shape = Unit)) +
  geom_line() +
  scale_y_log10(breaks = 10^seq(-12, 4)) +
  scale_shape_manual(values = 1:14) +
  annotation_logticks(sides = "l") +
  coord_cartesian(ylim = c(0.01, 5e4)) +
  facet_wrap(~Unit) +
  xlim(0, NA) +
  theme_bw()


out2 <- bind_rows(out_sep, out)

ggplot(out2, aes(uID, Conc, colour = Method)) +
  geom_pointrange(aes(ymin = LowHDI, ymax = TopHDI, shape = factor(pID)), position = position_dodge(0.4), size = 0.2) +
  scale_shape_manual(values = 1:9) +
  scale_y_log10()