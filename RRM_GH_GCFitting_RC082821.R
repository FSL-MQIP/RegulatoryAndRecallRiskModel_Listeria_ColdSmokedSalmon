## ------------------------------------------------------------------------
#  Recall risk model (for LM in CSS)
## ------------------------------------------------------------------------
## Script purpose: template codes for fitting and comparing different growth models.
## ------------------------------------------------------------------------
library(tidyverse);library(lmerTest);library(data.table);library(knitr);library(emmeans)
library(patternplot);library(car);library(multcomp);library(plyr);library(ggpubr);library(stringi)
library(nlsMicrobio);library(minpack.lm);library(AICcmodavg); library(grid)
## ------------------------------------------------------------------------
## Define primary growth models and die-off & regrowth models.

# Buchanan.
buchanan_log10 <- LOG10N ~ LOG10N0 + 
  (t >= lag) * (t <= (lag + (LOG10Nmax - LOG10N0) / mumax)) * mumax * (t - lag) + 
  (t >= lag) * (t > (lag + (LOG10Nmax - LOG10N0) / mumax)) * (LOG10Nmax - LOG10N0)

# Buchanan without lag.
buchanan_nl_log10 <- LOG10N ~ LOG10N0 + (t <= ((LOG10Nmax - LOG10N0) /mumax)) * mumax * t/ + (t > ((LOG10Nmax - LOG10N0) /mumax)) * (LOG10Nmax - LOG10N0)

# Gompertzm.
gompertzm_log10 <- LOG10N ~ LOG10N0 + (LOG10Nmax - LOG10N0) * exp(-exp(mumax * exp(1) * (lag - t)/((LOG10Nmax - LOG10N0)) + 1))

# Baranyi.
baranyi_log10 <- LOG10N ~ LOG10Nmax + log10((-1 + exp((mumax*log(10)) * lag) + exp((mumax*log(10)) * t))/(exp((mumax*log(10)) * t) - 1 + exp((mumax*log(10)) * lag) * 10^(LOG10Nmax - LOG10N0)))

# Baranyi without lag.
baranyi_nl_log10 <- LOG10N ~ (LOG10Nmax - log10(1 + (10^(LOG10Nmax - LOG10N0) - 1) * exp(-(mumax*log(10)) * t)))


# Weibull-Gompertzm.
WeiGom_tc <- LOG10N ~ LOG10N0 + 
  (t < tc) * (-(t/delta)^p) +
  (t >= tc) * ((-(tc/delta)^p)* (lag >= tc) + (LOG10Nmax - (LOG10N0 - ((tc/delta)^p)*(lag >= tc))) * exp(-exp(mumax * exp(1) * (lag - tc - t)/((LOG10Nmax - (LOG10N0 - ((tc/delta)^p)*(lag >= tc)))) + 1)))

# Weibull-Buchanan
WeiBuc_tc <- LOG10N ~ LOG10N0 + (t < tc) * (-(t/delta)^p) + 
  (t >= tc) * (t < lag) * (-(tc/delta)^p) + 
  (t >= lag) * (t <= (lag + (LOG10Nmax - LOG10N0 + (tc/delta)^p)/mumax)) * (-(tc/delta)^p * ((lag > tc) | (t > tc)) + mumax * (t - lag)) + 
  (t >= lag) * (t > (lag + (LOG10Nmax - LOG10N0 + (tc/delta)^p)/mumax)) * (LOG10Nmax - LOG10N0)

# Weibull-Buchanan without lag
WeiBuc_no_lag_tc <- LOG10N ~ LOG10N0 + (t < tc) * (-(t/delta)^p) + 
  (t >= tc) * (t <= (tc + (LOG10Nmax - LOG10N0 + (tc/delta)^p)/mumax)) * (-(tc/delta)^p + mumax * (t - tc)) + 
  (t > (tc + (LOG10Nmax - LOG10N0 + (tc/delta)^p)/mumax)) * (LOG10Nmax - LOG10N0)

# Weibull-Baranyi
WeiBar_tc <- LOG10N ~ 
  (t < tc) * (LOG10N0 - (t/delta)^p) +
  (t >= tc) * (lag >= tc) * (LOG10Nmax + log10((-1 + exp(mumax * (lag-tc)) + exp(mumax * t))/(exp(mumax * t) - 1 + exp(mumax * (lag-tc)) * 10^(LOG10Nmax - (LOG10N0 - (tc/delta)^p)))))


# Weibull-Baranyi without lag
WeiBar_no_lag_tc <- LOG10N ~ 
  (t < tc) * (LOG10N0 -(t/delta)^p) +
  (t >= tc) * ((LOG10Nmax - log10(1 + (10^(LOG10Nmax - LOG10N0 + ((tc/delta)^p)) - 1) * exp(-mumax * (t-tc)))))

## ------------------------------------------------------------------------
## Using the growth curve data of strain l3051 on cold-smoked salmon without nisin treatments as an example.

# Load in data.
l3051_C <- read.delim("l3051_C_withinoc_RC091619.txt")

# Five candidate growth models
# 1. Buchanan
# 2. Buchanan without lag
# 3. Gompertzm
# 4. Baranyi
# 5. Baranyi without lag

# Fit model 1: Buchanan.
# Preview the curve for setting up starting values of the parameters.
preview(buchanan_log10, l3051_C, list(LOG10N0 = 6,lag = 0.05, mumax = 0.20, LOG10Nmax = 8.5))
# Fit the data with the model.
fitl3051_C.buc_LM <- nlsLM(buchanan_log10, l3051_C, trace=T, 
                           list (LOG10N0 = 6,lag = 0.05, mumax = 0.20, LOG10Nmax = 8.5), 
                           control = nls.control(maxiter = 100, minFactor = 1/4096, warnOnly = T),
                           lower = c(LOG10N0 = 0,lag = 0, mumax = 0, LOG10Nmax = 0))
# Check the model fitting and the estimated values of the parameters.
plotfit(fitl3051_C.buc_LM)

# Fit model 2: Buchanan without lag.
preview(buchanan_nl_log10, l3051_C, list(LOG10N0 = 6,mumax = 0.2, LOG10Nmax = 8.5))
# Fit the data with the model.
fitl3051_C.buc_nl_LM <- nlsLM(buchanan_nl_log10, l3051_C, trace=T, 
                              list (LOG10N0 = 6,mumax = 0.2, LOG10Nmax = 8.5), 
                              control = nls.control(maxiter = 100, minFactor = 1/4096, warnOnly = T),
                              lower = c(LOG10N0 = 0,mumax = 0, LOG10Nmax = 0))
# Check the model fitting and the estimated values of the parameters.
plotfit(fitl3051_C.buc_nl_LM)

## Repeat this for the rest three models.

# Fit model 3: Gompertzm.
preview(gompertzm_log10, l3051_C, list(LOG10N0 = 6,lag = 0.05, mumax = 0.2, LOG10Nmax = 8.5))
fitl3051_C.gom_LM <- nlsLM(gompertzm_log10, l3051_C, trace=T, 
                           list (LOG10N0 = 6,lag = 0.05, mumax = 0.2, LOG10Nmax = 8.5), 
                           control = nls.control(maxiter = 100, minFactor = 1/4096, warnOnly = T),
                           lower = c(LOG10N0 = 0,lag = 0, mumax = 0, LOG10Nmax = 0))
plotfit(fitl3051_C.gom_LM)
# Fit model 4: Baranyi.
preview(baranyi_log10, l3051_C, list(LOG10N0 = 6,lag = 0.05, mumax = 0.2, LOG10Nmax = 8.5))
fitl3051_C.bar_LM <- nlsLM(baranyi_log10, l3051_C, trace=T, 
                           list (LOG10N0 = 6,lag = 0.1, mumax = 0.2, LOG10Nmax = 8.5), 
                           control = nls.control(maxiter = 100, minFactor = 1/4096, warnOnly = T),
                           lower = c(LOG10N0 = 0,lag = 0, mumax = 0, LOG10Nmax = 0))
plotfit(fitl3051_C.bar_LM)
# Fit model 5: Baranyi without lag.
preview(baranyi_nl_log10, l3051_C, list(LOG10N0 = 6, mumax = 0.2, LOG10Nmax = 8.5))
fitl3051_C.bar_nl_LM <- nlsLM(baranyi_nl_log10, l3051_C, trace=T, 
                              list (LOG10N0 = 6, mumax = 0.2, LOG10Nmax = 8.5), 
                              control = nls.control(maxiter = 100, minFactor = 1/4096, warnOnly = T),
                              lower = c(LOG10N0 = 0, mumax = 0, LOG10Nmax = 0))
plotfit(fitl3051_C.bar_nl_LM)

# Compare across models to select the best-fit model.
candidate_models <- list()
candidate_models[[1]] <- fitl3051_C.buc_LM
candidate_models[[2]] <- fitl3051_C.buc_nl_LM
candidate_models[[3]] <- fitl3051_C.gom_LM
candidate_models[[4]] <- fitl3051_C.bar_LM
candidate_models[[5]] <- fitl3051_C.bar_nl_LM
mod.names <- c("Buchanan","Buchanan_withou_lag", "Gompertzm", "Baranyi", "Baranyi_without_lag")
title_string <- paste("Isolate ", "l3-051", "0 ppm", sep=" ")
title_string
output_bic <- bictab(cand.set = candidate_models, modnames = mod.names, sort = TRUE)
print(title_string)
print(output_bic)
## ------------------------------------------------------------------------




