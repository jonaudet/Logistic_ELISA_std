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

unkn <- read_csv("Data.csv") %>%
  filter(pID == 3) %>%
  mutate(uID = as.numeric(factor(Samp)),
         pID = pID - 2) %>%
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


initial <- function(N, N_plates, N_grp){
  inits <- list(pred_std_raw = rnorm(1, 0, 1),
                sigma = abs(rnorm(1, 0, 1)),
                mu_Bottom = abs(rnorm(1, 0.05, 0.02)),
                mu_Span = rnorm(1, 3.5, 0.1),
                mu_log_Inflec = rnorm(1, 0, 1),
                mu_Slope = abs(rnorm(1, 1, 0.5)),
                log_theta = runif(N_grp - 1, -5, 6),
                sigma_x = rexp(1, 1),
                log_x_raw = rnorm(N, 0, 1))
  return(inits)
}

# Run the model

inits <- lapply(1:8, function(x) initial(96, 1, max(unkn$uID)))

chains <- 1:8
xs <- exp(seq(-7, 4, length.out = 40))
dat <- expand.grid(chains, xs, KEEP.OUT.ATTRS = F)
inv_logit <- function(x) return(1 / (1 + exp(-x)))
bottom <- sapply(1:8, function(i) inits[[i]]$mu_Bottom)
span <- sapply(1:8, function(i) inits[[i]]$mu_Span)
log_Inf <- sapply(1:8, function(i) inits[[i]]$mu_log_Inflec)
slope <- sapply(1:8, function(i) inits[[i]]$mu_Slope)
sigma <- sapply(1:8, function(i) inits[[i]]$sigma)
sigma_x <- sapply(1:8, function(i) inits[[i]]$sigma_x)
sigma_ratio <- sigma/sigma_x

sigma_ratio

dat <- mutate(dat,
              OD = bottom[Var1] + span[Var1] * inv_logit((log(Var2) - log_Inf[Var1]) * slope[Var1]),
              Var1 = factor(Var1)) %>%
  select(Chain = Var1, Conc = Var2, OD)

mutate(unkn, Conc = 4500 * Dilution) %>%
  filter(Std == T) %>%
  ggplot(aes(x = Conc, y = OD, colour = Std, group = uID)) +
  geom_point(alpha = 0.4) +
  geom_line(data = dat, aes(colour = Chain, group = Chain)) +
  scale_x_log10() +
  theme_bw()

res3 <- stan(file = "logistic_OD_4p_UnknOnly.stan",
             data = list(N_unkn = nrow(unkn),
                         N_unkn_grp = max(unkn$uID),
                         uID = unkn$uID,
                         Unknown = unkn$OD,
                         ser_dils = ser_dilutions,
                         mu_Std = 4500,
                         sigma_std = 200),
             init = inits, chains = 8,
             iter = 14000, warmup = 10000, refresh = 200, control = list(adapt_delta = 0.95))#, max_treedepth = 15))

stan_ess(res3)
stan_ac(res3)
stan_rhat(res2)

traceplot(res3, pars = "mu_Span", inc_warmup = T) + facet_wrap(~chain, ncol = 1)

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
