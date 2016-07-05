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

dat <- LoadData(Folder = "E:/Jonathan/Bello_GP_IgG", DilFile = "Dilution.xlsx", DataFile = "Data.xlsx", LayoutFile = "Layout.xlsx") %>%
  separate(Plate, c("Noun", "pID"), convert = TRUE) %>%
  mutate(Plate = paste(Noun, pID)) %>%
  filter(!is.na(OD)) %>%
  select(-Noun) %>%
  filter(pID >= 12) %>%
  mutate(pID = pID - 11)

unkn <- dat %>%
  #filter(Samp != "std") %>%
  mutate(uID = as.numeric(factor(Samp))) %>%
  arrange(uID, Dilution) %>%
  mutate(Dilution = 1 / Dilution)

dat <- filter(dat, Samp == "std") %>%
  arrange(pID, Dilution)

start_dils <- unkn %>%
  mutate(uID = factor(uID)) %>%
  group_by(uID) %>%
  dplyr::summarise(SDil = max(Dilution)) %>%
  ungroup %>%
  .$SDil

ser_dilutions <- unkn %>%
  .$Dilution

std_ser_dilutions <- dat %>%
  mutate(Dilution = 1 / Dilution) %>%
  .$Dilution


res2 <- stan(file = "logistic_OD_5p_multiUnkn.stan",
             data = list(N_unkn = nrow(unkn),
                         N_unkn_grp = max(unkn$uID),
                         uID = unkn$uID,
                         N_plates = 7,
                         Unknown = unkn$OD,
                         pID = unkn$pID,
                         ser_dils = ser_dilutions,
                         mu_Std = 4500,
                         sigma_std = 200),
             iter = 8000, warmup = 4000, refresh = 200, control = list(adapt_delta = 0.95)) #iter = 15000, warmup = 5000, thin = 10, control = list(adapt_delta = 0.95))

ggplot(dat, aes(x = 4500 * Dilution, y = OD, colour = Plate)) +
  geom_point() +
  scale_x_log10() +
  facet_wrap(~Plate, ncol = 3)

print(res2, pars = c("cOD",
                     "unkn_cOD",
                     "sigma_indiv",
                     "unkn_sigma",
                     "log_x_init",
                     "x_init"), include = F)

dat2 <- dat

cOD <- rstan::extract(res2, "unkn_cOD")$unkn_cOD
dat2$cOD <- apply(cOD, 2, median)
dat2 <- bind_cols(dat2,
                 ldply(apply(cOD, 2, HDI),
                       function(x) return(x)))
dat2$Conc <- 4500/dat2$Dilution

output <- unkn

cOD <- rstan::extract(res2, "x_init")$x_init
output$Median <- apply(cOD, 2, median)
errors <- ldply(apply(cOD, 2, HDI),
                function(x) return(x))
output$TopHDI <- errors$HigherHDI
output$LowHDI <- errors$LowerHDI
output$Conc <- output$Median


ggplot(dat2, aes(Conc, OD)) +
  # geom_point(alpha = 0.4) +
  #   geom_line(colour = "blue") +
  #   geom_line(aes(y = cOD)) +
  scale_x_log10(breaks = 10^seq(floor(log10(min(output$Conc))), ceiling(log10(max(output$Conc))), by = 1)) +
  coord_cartesian(xlim = c(1e-4, 35)) +
  # facet_wrap(~Plate, ncol = 2) +
  # geom_pointrange(aes(y = cOD, ymin = LowHDI, ymax = TopHDI), colour = "red", size = 0.3) +
  geom_point(data = output, aes(Conc, OD, colour = factor(uID))) +
  geom_errorbarh(data = output, aes(y = OD, xmin = LowHDI, xmax = TopHDI, colour = factor(uID))) +
  scale_colour_discrete(guide = "none") +
  facet_wrap(~pID, ncol = 3)

out <- unkn %>%
  filter(Samp != "std") %>%
  arrange(uID) %>%
  group_by(uID) %>%
  top_n(1, OD) %>%
  ungroup %>%
  separate(Samp, c("Group", "Samp"), sep = "-") %>%
  separate(Samp, c("Animal", "Day"), sep = "_") %>%
  mutate(Day = as.numeric(Day))

theta <- rstan::extract(res2, "theta")$theta
out$Conc <- apply(theta, 2, median)
errors <- ldply(apply(theta, 2, HDI),
                function(x) return(x))
out$TopHDI <- errors$HigherHDI
out$LowHDI <- errors$LowerHDI

ggplot(out, aes(Day, Conc, colour = Group, fill = Animal)) +
  geom_pointrange(aes(ymin = LowHDI, ymax = TopHDI, shape = Animal)) +
  geom_line() +
  scale_y_log10(breaks = 10^seq(-5, 4)) +
  annotation_logticks(sides = "l") +
  coord_cartesian(ylim = c(0.01, 5e3)) +
  xlim(0, NA) +
  theme_bw()


sessionInfo()

ldply(apply(rstan::extract(res2, "theta")$theta, 2, HDI), function(x) return(x))
