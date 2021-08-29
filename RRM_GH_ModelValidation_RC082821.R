## ------------------------------------------------------------------------
#  Recall risk model (for LM in CSS)
## ------------------------------------------------------------------------
## Script purpose: Model validation
## ------------------------------------------------------------------------
library(tidyverse);library(truncnorm);library(Dict);library(Matching);library(NCmisc);library(quanteda);library(ggpubr)
## ------------------------------------------------------------------------
cat(strftime(Sys.time(),"%Y-%m-%d %H:%M:%S", tz = ""),'\n'," package import",'\n')
## ------------------------------------------------------------------------
## Define primary growth models and die-off & regrowth models.
# Model 1: Buchanan without lag (in LOG10 scale)
buchanan_nl_log10 = function(t,mumax,LOG10N0,LOG10Nmax){
  LOG10N <- LOG10N0 + (t <= ((LOG10Nmax - LOG10N0) /mumax)) * mumax * t + (t > ((LOG10Nmax - LOG10N0) /mumax)) * (LOG10Nmax - LOG10N0)
  return(LOG10N)
}
# Model 2: Baranyi without lag (in LOG10 scale)
baranyi_nl_log10 = function(t,mumax,LOG10N0,LOG10Nmax){
  LOG10N <- (LOG10Nmax - log10(1 + (10^(LOG10Nmax - LOG10N0) - 1) * exp(-(mumax*log(10)) * t)))
  return(LOG10N)
}
# Model 3: Weibull-Gompertzm (in LOG10 scale)
WeiGom_tc_II = function(t,p,delta,lag,mumax,LOG10N0,LOG10Nmax,tc) {
  LOG10N <- LOG10N0 + 
    (t < tc) * (-(t/delta)^p) +
    (t >= tc) * 
    ((-(tc/delta)^p)* (lag >= tc) + (LOG10Nmax - (LOG10N0 - ((tc/delta)^p)*(lag >= tc))) * exp(-exp(mumax * exp(1) * (lag - tc - t)/((LOG10Nmax - (LOG10N0 - ((tc/delta)^p)*(lag >= tc)))) + 1)))
  return(LOG10N)
}
# Model 4: Weibull-Buchanan without lag (in LOG10 scale)
WeiBuc_nl_tc = function(t,p,delta,mumax,LOG10N0,LOG10Nmax,tc){
  LOG10N <- LOG10N0 + 
    (t < tc) * (-(t/delta)^p) + 
    (t >= tc) * (t <= (tc + (LOG10Nmax - LOG10N0 + (tc/delta)^p)/mumax)) * (-(tc/delta)^p + mumax * (t - tc)) + 
    (t > (tc + (LOG10Nmax - LOG10N0 + (tc/delta)^p)/mumax)) * (LOG10Nmax - LOG10N0)
  return(LOG10N)
}
# Model 5: Weibull-Baranyi (in LOG10 scale)
WeiBar_tc_II = function(t,p,delta,lag,mumax,LOG10N0,LOG10Nmax,tc) {
  LOG10N <- (t < tc) * (LOG10N0 - (t/delta)^p) + (t >= tc) * (lag >= tc) * (LOG10Nmax + log10((-1 + exp((mumax*log(10)) * (lag-tc)) + exp((mumax*log(10)) * t))/(exp((mumax*log(10)) * t) - 1 + exp((mumax*log(10)) * (lag-tc)) * 10^(LOG10Nmax - (LOG10N0 - (tc/delta)^p)))))
  return(LOG10N)
}
# Model 6: Weibull-Baranyi without lag (in LOG10 scale)
WeiBar_nl_tc = function(t,p,delta,mumax,LOG10N0,LOG10Nmax,tc){
  LOG10N <- (t < tc) * (LOG10N0 - (t/delta)^p) + (t >= tc) * ((LOG10Nmax -  log10(1 + (10^(LOG10Nmax - LOG10N0 + ((tc/delta)^p)) - 1) * exp(-(mumax*log(10)) * (t - tc)))))
  return(LOG10N)
}
# Model 7: Buchanan with lag (in log10 scale)
buchanan_log10 = function(t,mumax,LOG10N0,LOG10Nmax,lag){
  LOG10N <- LOG10N0 + #base population
    (t >= lag) * (t <= (lag + (LOG10Nmax - LOG10N0) / mumax)) * #if in growth phase
    mumax * (t - lag) + #log-linear growth (assume positive)
    (t >= lag) * (t > (lag + (LOG10Nmax - LOG10N0) / mumax)) * # if in stationary phase (first condition unnecessary)
    (LOG10Nmax - LOG10N0) #take value of Nmax
  return(LOG10N)
}
# Model 8: Baranyi with lag (in log10 scale)
# (mumax*log(10))
baranyi_log10 = function(t,mumax,LOG10N0,LOG10Nmax,lag){
  LOG10N <- LOG10Nmax + log10((-1 + exp((mumax*log(10)) * lag) + exp((mumax*log(10)) * t))/(exp((mumax*log(10)) * t) - 1 + exp((mumax*log(10)) * lag) * 10^(LOG10Nmax - LOG10N0)))
  return(LOG10N)
}

# Generic function for implementing primary growth models and die-off & regrowth models.
LOG10N_calculation <- function(model_name = "buchanan", t,p,delta,lag,mumax,LOG10N0,LOG10Nmax,tc) {
  if (model_name == "buchanan_without_lag") {
    return(buchanan_nl_log10(t,mumax,LOG10N0,LOG10Nmax))
  } else if (model_name == "baranyi_without_lag") {
    return(baranyi_nl_log10(t,mumax,LOG10N0,LOG10Nmax))
  } else if (model_name == "weibull_buchanan_without_lag") {
    return(WeiBuc_nl_tc(t,p,delta,mumax,LOG10N0,LOG10Nmax,tc))
  } else if (model_name == "weibull_baranyi") {
    return(WeiBar_tc_II(t,p,delta,lag,mumax,LOG10N0,LOG10Nmax,tc))
  } else if (model_name == "weibull_baranyi_without_lag") {
    return(WeiBar_nl_tc(t,p,delta,mumax,LOG10N0,LOG10Nmax,tc))
  } else if (model_name == "buchanan") {
    return(buchanan_log10(t,mumax,LOG10N0,LOG10Nmax,lag))
  } else if (model_name == "baranyi") {
    return(baranyi_log10(t,mumax,LOG10N0,LOG10Nmax,lag))
  } else {
    stop(paste0(model_name, " is not a valid model name. Must be one of buchanan, baranyi, gompertz"))
  }
}
## ------------------------------------------------------------------------
## Define other equations for use in the model.
# Function for converting Mumax from reference to actual temperature.
muAtNewTemp_CSS <- function(newTemp, oldMu, oldTemp = 5, T0 = -2.86) {
  numerator <- newTemp - T0
  denom <- oldTemp - T0
  newMu <- ((numerator / denom)^2) * oldMu
  return(newMu)
}

# Ramdom sampling from truncated weibull distribution.
rtweibull <- function(n, shape, scale, a=0, b=Inf)
{
  stopifnot(n > 0 & all(scale > 0) & all(shape > 0))
  x <- runif(n)
  Fa <- pweibull(a, shape, scale)
  Fb <- pweibull(b, shape, scale)
  y <- (1-x)*Fa + x*Fb
  return(qweibull(y, shape, scale))
}

# Getting quantiles of truncated weibull distributions.
qtweibull <- function(p, shape, scale, a=0, b=Inf)
{
  stopifnot( all(p >= 0 & p <= 1) & all( scale > 0 ) & all(shape > 0) )
  Fa <- pweibull(a, shape, scale)
  Fb <- pweibull(b, shape, scale)
  pNew <- p * (Fb - Fa) + Fa
  x <- qweibull(pNew, shape, scale)
  return(x)
}

cat(strftime(Sys.time(),"%Y-%m-%d %H:%M:%S", tz = ""),'\n',"functions defined",'\n')
## ------------------------------------------------------------------------
## Get experimental data.
HighInocData <- read.csv("ValidationData_HighInoc_RC030421.csv")
HI_N0D15 <- HighInocData %>% filter(Nisin == "Minus", Day == 15, Condition == "BHI")
HI_N0D30 <- HighInocData %>% filter(Nisin == "Minus", Day == 30, Condition == "BHI")
HI_N25D15 <- HighInocData %>% filter(Nisin == "Plus", Day == 15, Condition == "BHI")
HI_N25D30 <- HighInocData %>% filter(Nisin == "Plus", Day == 30, Condition == "BHI")
## ------------------------------------------------------------------------
## Set up model-associated parameter values and distributions.

# Growth parameters.
GrowthPars <- dict("0" = read.delim("growth_parameter_0ppm"), "5" = read.delim("growth_parameter_5ppm"),
                   "10" = read.delim("growth_parameter_10ppm"), "20" = read.delim("growth_parameter_20ppm"))
# Strains.serotype_vec <- c("12a", "12b", "4b")
strain_dict <- dict("0" = c("c1-111", "f2-237", "f6-366", "l3-051"), "5" = c("c1-111", "f6-366", "l3-051"),
                    "10" = c("c1-111", "f2-237", "f6-366", "l3-051"), "20" = c("c1-111", "f2-237", "f6-366", "l3-051"))
# Maximum population density.
LOG10Nmax <- dict("0" = dict(shape=41.341902, scale=9.126792), "5" = dict(shape=21.59984, scale=8.712891),
                  "10" = dict(shape=14.54644, scale=8.497603), "20" = dict(shape=8.438854, scale=8.208065))
# Minimum growth temperature.
mTmin <- rnorm(lot_iter_number, mean = -2.86, sd = 0.459)
ln_sdTmin <- rnorm(lot_iter_number, mean = 0.638, sd = 0.208); sdTmin <- exp(ln_sdTmin)
# Initial contamination level.
mLOG10N0 <- 5.9818182; sdLOG10N0 <- 0.1140429
# Maximum growth rate.
mmumax <-  rnorm(lot_iter_number, mean = 6.7091025, sd= 0.2596307)
ln_sdmumax <- rgamma(lot_iter_number, shape = 159.8772, scale = 0.005836614); sdmumax <- exp(ln_sdmumax)
# Random error of the growth models and the Die-off & Regrowth models.
mlnsigma <- -1.20
sdlnsigma <- 0.0185
## ---------------------------------------------------------------------------------------------------
## Model validation: 0 ppm nisin, day 15.

# Set basic model parameters.
lot_iter_number=10;package_number=10;
nisin_concentration=0; temp = 7; day=15
growth_parameter <- GrowthPars$get(as.character(nisin_concentration))

# Run the model.
source("RRM_GH_ModelValidation_MainModel_RC082821.R")

# Retrieve simulated data.
HIN0D15_SimulateDf <- simulated_df; HIN0D15_SimulateData <- simulated_data
# Retrieve experimental data.
HIN0D15_ExperimentalData <- HI_N0D15$Log_CFU

# Mann-Whitney U test.
MW_U_HIN0D15 <- wilcox.test(HIN0D15_ExperimentalData,HIN0D15_SimulateData, paired = FALSE)

# Calculate coverage rate.
Simulated_In_Experimental <- ifelse(HIN0D15_SimulateData >= 8.17 & HIN0D15_SimulateData <= 9.27, 1, 0)
Coverage_Rate <- sum(Simulated_In_Experimental)/length(Simulated_In_Experimental)

## ---------------------------------------------------------------------------------------------------
## Model validation: 25 ppm nisin, day 15.

# Set basic model parameters.
lot_iter_number=10;package_number=10;
nisin_concentration=20; temp = 7; day=15
growth_parameter <- GrowthPars$get(as.character(nisin_concentration))

# Run the model.
source("ModelValidation_MainModel_RC030421.R")
# Retrieve simulated data.
HIN20D15_SimulateDf <- simulated_df; HIN20D15_SimulateData <- simulated_data
# Retrieve experimental data.
HIN25D15_ExperimentalData <- HI_N25D15$Log_CFU

# Mann-Whitney U test.
MW_U_HIN20D15 <- wilcox.test(HIN25D15_ExperimentalData,HIN20D15_SimulateData, paired = FALSE)
# Calculate coverage rate.
Simulated_In_Experimental <- ifelse(HIN20D15_SimulateData >= 4.88 & HIN20D15_SimulateData <= 8.52, 1, 0)
Coverage_Rate <- sum(Simulated_In_Experimental)/length(Simulated_In_Experimental)
## ------------------------------------------------------------------------


