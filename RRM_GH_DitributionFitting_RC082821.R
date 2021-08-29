## ------------------------------------------------------------------------
#  Recall risk model (for LM in CSS)
## ------------------------------------------------------------------------
## Script purpose: template codes for fitting and comparing parametric distributions.
## ------------------------------------------------------------------------
library(tidyverse);library(lmerTest);library(data.table);library(knitr);library(emmeans)
library(patternplot);library(car);library(multcomp);library(plyr);library(ggpubr);library(stringi)
library(nlsMicrobio);library(minpack.lm);library(AICcmodavg)
library(fitdistrplus); library(mc2d)
## ------------------------------------------------------------------------
## Using initial contamination level as an example.
## ------------------------------------------------------------------------
# Load in data.
init_contam_df <- read.csv('initial_contam_afssa.csv')
## ------------------------------------------------------------------------
## Determine the distribution of N0.

# Check the cumulative distribution of the censored data.
plotdistcens(init_contam_df, NPMLE = FALSE)

# Fit candidate distributions: Weibull, Gamma, and Lognormal.
# Weibull.
fd_init_contam_weibull <- fitdistcens(init_contam_df, "weibull")
plot(fd_init_contam_weibull)
# Gamma.
fd_init_contam_gamma <- fitdistcens(init_contam_df, "gamma")
plot(fd_init_contam_gamma)
# Lognormal.
fd_init_contam_lnorm <- fitdistcens(init_contam_df, "lnorm")
plot(fd_init_contam_lnorm)

# Plot the fitting of different distributions to identify the best fit.
cdfcompcens(list(fd_init_contam_weibull, fd_init_contam_gamma, fd_init_contam_lnorm), 
            legendtext = c("weibull",  "gamma", "lognormal"))
## ------------------------------------------------------------------------
## Use the bootstrap resampling method to determine the distributions of mLnN0 and sdLnN0.
boot.init_contam <- bootdistcens(fd_init_contam_lnorm, niter = 1001)
plot(boot.init_contam)

# Fit a distribution for mLnN0.
init_contam_meanlog <- boot.init_contam$estim$meanlog
plotdist(init_contam_meanlog, histo = TRUE, demp = TRUE)
descdist(init_contam_meanlog, boot = 1000)
qqp(init_contam_meanlog, "norm")
# In this case, normal distribution should be the best fit.
# Fit a normal distribution
init_contam_meanlog_norm <- fitdist(init_contam_meanlog, "norm")
summary(init_contam_meanlog_norm)


# Fit a distribution for sdLnN0.
init_contam_sdlog <- boot.init_contam$estim$sdlog
plotdist(init_contam_sdlog, histo = TRUE, demp = TRUE)
descdist(init_contam_sdlog, boot = 1000)
qqp(init_contam_sdlog, "norm")
# In this case, normal distribution should be the best fit.
# Fit a normal distribution.
init_contam_sdlog_norm <- fitdist(init_contam_sdlog, "norm")
summary(init_contam_sdlog_norm)
## ------------------------------------------------------------------------


