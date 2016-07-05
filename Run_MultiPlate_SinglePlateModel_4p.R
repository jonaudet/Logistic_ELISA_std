if(!require(checkpoint)){
  install.packages("checkpoint")
  library(checkpoint)
}
checkpoint("2016-01-13")
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
    dplyr::mutate(Well = gsub("_", "", Well))

  return(Data)
}

unkn <- read_csv("Data.csv") %>%
  filter(Plate != "Plate 12") %>%
  mutate(uID = as.numeric(factor(Samp)),
         pID = as.numeric(factor(Plate))) %>%
  arrange(uID, Dilution) %>%
  mutate(Std = Samp == "std")

ser_dilutions <- unkn %>%
  .$Dilution

# Plot of each plate's standard and corresponding unknowns
mutate(unkn, Conc = 4500 * Dilution) %>%
  ggplot(aes(x = Conc, y = OD, colour = Std, group = uID)) +
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


initial <- function(N, N_plates, N_grp){
  inits <- list(pred_std_raw = rnorm(1, 0, 1),
                sigma = abs(rnorm(1, 0, 1)),
                mu_Bottom = abs(rnorm(N_plates, 0.05, 0.02)),
                mu_Span = rnorm(N_plates, 3.5, 0.1),
                mu_log_Inflec = rnorm(N_plates, 0, 1),
                mu_Slope = abs(rnorm(N_plates, 1, 0.5)),
                log_theta = runif(N_grp - 1, -5, 6),
                sigma_x = rexp(1, 1),
                log_x_raw = rnorm(N, 0, 1))
  return(inits)
}

# Run the model

sep <- lapply(1:6, function(i){
  print(i)
  df <- filter(unkn, pID == i) %>%
    mutate(uID = as.numeric(factor(Samp)),
           pID = as.numeric(factor(Plate))) %>%
    arrange(uID, Dilution) %>%
    mutate(Std = Samp == "std")
  inits <- lapply(1:4, function(x) initial(nrow(df), 1, max(df$uID)))
  ser_dilutions <- df %>%
    .$Dilution
  res <- stan(file = "logistic_OD_4p_UnknOnly.stan",
               data = list(N_unkn = nrow(df),
                           N_unkn_grp = max(df$uID),
                           uID = df$uID,
                           Unknown = df$OD,
                           ser_dilutions = ser_dilutions,
                           mu_Std = 4500,
                           sigma_std = 200),
               init = inits, chains = 4,
               iter = 12000, warmup = 8000, refresh = 200, control = list(adapt_delta = 0.95))
  return(res)
})

output <- unkn
out_sep <- ldply(1:6, function(i){
  print(i)
  df <- filter(output, pID == i) %>%
    mutate(uID = as.numeric(factor(Samp)),
           pID = as.numeric(factor(Plate))) %>%
    arrange(uID, Dilution) %>%
    mutate(Std = Samp == "std")
  inits <- lapply(1:4, function(x) initial(nrow(df), 1, max(df$uID)))
  ser_dilutions <- df %>%
    .$Dilution

  cOD <- rstan::extract(sep[[i]], "x")$x
  df$Median <- apply(cOD, 2, median)
  errors <- ldply(apply(cOD, 2, HDI),
                  function(x) return(x))
  df$TopHDI <- errors$HigherHDI
  df$LowHDI <- errors$LowerHDI
  df$Conc <- df$Median
  return(df)
})

ggplot(out_sep, aes(Conc, OD)) +
  scale_x_log10(breaks = 10^seq(floor(log10(min(out_sep$Conc))), ceiling(log10(max(out_sep$Conc))), by = 1)) +
  coord_cartesian(xlim = c(1e-4, 35), ylim = c(0, 4)) +
  #geom_point(aes(colour = factor(uID))) +
  geom_text(aes(label = Samp, colour = factor(uID))) +
  geom_errorbarh(aes(xmin = LowHDI, xmax = TopHDI, colour = factor(uID))) +
  scale_colour_discrete(guide = "none") +
  facet_wrap(~Plate, ncol = 3)


out_sep <- ldply(1:6, function(i){
  print(i)
  out <- unkn %>%
  filter(Samp != "std", pID == i) %>%
    mutate(uID = as.numeric(factor(Samp)),
           pID = as.numeric(factor(Plate))) %>%
    arrange(uID) %>%
  group_by(uID) %>%
  top_n(1, OD) %>%
  ungroup %>%
  separate(Samp, c("Group", "Samp"), sep = "-") %>%
  separate(Samp, c("Unit", "Week"), sep = "_") %>%
  mutate(Week = as.numeric(Week))

  theta <- rstan::extract(sep[[i]], "theta")$theta
  out$Conc <- apply(theta, 2, median)
  errors <- ldply(apply(theta, 2, HDI),
                  function(x) return(x))
  out$TopHDI <- errors$HigherHDI
  out$LowHDI <- errors$LowerHDI
  return(out)
})

ggplot(out_sep, aes(Week, Conc, colour = Group, fill = Unit)) +
  geom_pointrange(aes(ymin = LowHDI, ymax = TopHDI, shape = Unit)) +
  geom_line() +
  scale_y_log10(breaks = 10^seq(-12, 4)) +
  annotation_logticks(sides = "l") +
  coord_cartesian(ylim = c(0.1, 5e3)) +
  xlim(0, NA) +
  theme_bw()

ser_dilutions <- unkn %>%
  .$Dilution

inits <- lapply(1:4, function(x) initial(nrow(unkn), max(unkn$pID), max(unkn$uID)))

res2 <- stan(file = "logistic_OD_4p_MultiPlate.stan",
             data = list(N_unkn = nrow(unkn),
                         N_unkn_grp = max(unkn$uID),
                         uID = unkn$uID,
                         Unknown = unkn$OD,
                         ser_dils = ser_dilutions,
                         mu_Std = 4500,
                         sigma_std = 200,
                         N_plates = max(unkn$pID),
                         pID = unkn$pID),
             init = inits, chains = 4,
             iter = 14000, warmup = 10000, refresh = 200, control = list(adapt_delta = 0.95))#, max_treedepth = 15))

stan_ess(res2)
stan_ac(res2)
stan_rhat(res2)

# Look at the curves and the unknowns, the standards are shown in red all other colors are unknowns

output <- unkn

cOD <- rstan::extract(res2, "x")$x
output$Median <- apply(cOD, 2, median)
errors <- ldply(apply(cOD, 2, HDI),
                function(x) return(x))
output$TopHDI <- errors$HigherHDI
output$LowHDI <- errors$LowerHDI
output$Conc <- output$Median

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
  separate(Samp, c("Group", "Samp"), sep = "-") %>%
  separate(Samp, c("Unit", "Day"), sep = "_") %>%
  mutate(Day = as.numeric(Day))

theta <- rstan::extract(res2, "theta")$theta
out$Conc <- apply(theta, 2, median)
errors <- ldply(apply(theta, 2, HDI),
                function(x) return(x))
out$TopHDI <- errors$HigherHDI
out$LowHDI <- errors$LowerHDI

ggplot(out, aes(Day, Conc, colour = Group, fill = Unit)) +
  geom_pointrange(aes(ymin = LowHDI, ymax = TopHDI, shape = Unit)) +
  geom_line() +
  scale_y_log10(breaks = 10^seq(-12, 4)) +
  annotation_logticks(sides = "l") +
  coord_cartesian(ylim = c(0.00001, 5e3)) +
  xlim(0, NA) +
  theme_bw()
