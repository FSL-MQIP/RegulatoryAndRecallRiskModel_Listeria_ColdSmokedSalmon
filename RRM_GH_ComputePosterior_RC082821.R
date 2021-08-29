## ------------------------------------------------------------------------
#  Recall risk model (for LM in CSS)
## ------------------------------------------------------------------------
## Script purpose: update the prior distributions of mMumax and sdMumax with experimental data.
## ------------------------------------------------------------------------
library(NCmisc);library(BayesianTools);library(tidyverse);library(fitdistrplus)
## ------------------------------------------------------------------------
## Define functions to use.
# Convert mumax from the reference to actual temperature.
muAtNewTemp_CSS <- function(newTemp, oldMu, oldTemp = 5, T0 = -2.86) {
  numerator <- newTemp - T0
  denom <- oldTemp - T0
  newMu <- ((numerator / denom)^2) * oldMu
  return(newMu)
}
## ------------------------------------------------------------------------

##-------------------------------------------------------
# Prior distributions (Delignette-Muller et al., 2006)
# mmumax ~ gamma(shape=69.7, scale=0.0896)
# log(sdmumax) ~ norm(-1.03, 0.191)
##-------------------------------------------------------
## Read in new Mumax data.
Ji_mumax <- read.csv("mumax_Jihun_RC022821.csv")
Silin_mumax <- read.csv("mumax_Silin_RC022821.csv")
LI_mumax <- read.csv("mumax_LI_RC022821.csv")
##-------------------------------------------------------
## Identify outliers in the data sets.
# Jihun's data.
Outliers_Ji <- boxplot.stats(Ji_mumax$mumax_25C)$out
OutliersRows_Ji <- which(Ji_mumax$mumax_25C %in% c(Outliers_Ji))
Ji_mumax_Filtered <- Ji_mumax[-OutliersRows_Ji, ]
# Silin's data.
Outliers_Silin <- boxplot.stats(Silin_mumax$mumax_25C)$out
OutliersRows_Silin <- which(Silin_mumax$mumax_25C %in% c(Outliers_Silin))
Silin_mumax_Filtered <- Silin_mumax[-OutliersRows_Silin, ]
# LI data.
Outliers_LI <- boxplot.stats(LI_mumax$LI_mumaxLN25)$out
LI_mumax_Filtered <- LI_mumax
##-------------------------------------------------------
## Combine all of the new data.
y = c(Ji_mumax_Filtered$mumax_25C, Silin_mumax_Filtered$mumax_25C, LI_mumax_Filtered$LI_mumaxLN25)
##-------------------------------------------------------
## Set up functions for running MCMC.
likelihood <- function(param) {
  a = param[1]
  b = param[2]
  singlelikelihoods = dnorm(y, mean = a, sd = exp(b), log = T)
  sumll = sum(singlelikelihoods)
  return(sumll)
}

density <- function(param){
  a = param[1]
  b = param[2]
  adensity = dgamma(a, shape=69.7, scale=0.0896, log = T)
  bdensity = dnorm(b, mean=1.03, sd = 0.191, log = T)
  return(adensity+bdensity)
}

sampler = function(n=1){
  asampler = rgamma(n, shape=69.7, scale=0.0896)
  bsampler = rnorm(n, mean=1.03, sd = 0.191)
  return(cbind(asampler,bsampler))
}

posterior <- function(param){
  return (likelihood(param) + prior(param))
}

prior <- createPrior(density = density, sampler = sampler,
                     lower = c(-10,-10), upper = c(10,10), best = NULL)
##-------------------------------------------------------
## Run MCMC with the "Metropolis" method.
set.seed(03042021)
bayesianSetup <- createBayesianSetup(likelihood = likelihood, prior = prior)
settings = list(iterations = 10000, nrChains = 3)
out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings, sampler = "Metropolis")
# Combine the results of three chains.
outChain1_df <- as.data.frame(out[[1]]$chain)
colnames(outChain1_df)[1:2] <- c("mmumax", "sdmumax")
outChain2_df <- as.data.frame(out[[2]]$chain)
colnames(outChain2_df)[1:2] <- c("mmumax", "sdmumax")
outChain3_df <- as.data.frame(out[[3]]$chain)
colnames(outChain3_df)[1:2] <- c("mmumax", "sdmumax")
Output_Combined <- rbind(outChain1_df[5001:10000,], outChain2_df[5001:10000,], outChain3_df[5001:10000,])
# Write out the results.
write.csv(Output_Combined, file = "PosteriorMCMCResults_mumax_RC031021V2.csv")
##-------------------------------------------------------
## Fit the MCMC results of mMumax and sdMumx with parametric distributions.
##-------------------------------------------------------
# Import files.
Posterior_MCMCResults <- read.csv("PosteriorMCMCResults_mumax_RC031021V2.csv")
mmumax_MCMCResults <- Posterior_MCMCResults$mmumax
sdmumax_MCMCResults <- Posterior_MCMCResults$sdmumax
##-------------------------------------------------------
# mMumax.
par(mar=c(2.0,2.0,2.0,2.0))
plotdist(mmumax_MCMCResults, histo = TRUE, demp = TRUE)
descdist(mmumax_MCMCResults, boot = 1000)

# Fit candidate parametric distributions: Weibull, Gamma, Lognormal, and Normal.
mmumax_weibull <- fitdist(mmumax_MCMCResults, "weibull")
mmumax_gamma <- fitdist(mmumax_MCMCResults, "gamma")
mmumax_lnorm <- fitdist(mmumax_MCMCResults, "lnorm")
mmumax_norm <- fitdist(mmumax_MCMCResults, "norm")

# Comparing different distributions - goodness-of-fit plots.
par(mfrow = c(2, 2))
par(mar=c(2.0,2.0,2.0,1.0))
plot.legend <- c("Weibull", "gamma", "lognormal", "norm")
# Red---Weibull; Green---gamma; Blue---lognormal; Light blue---normal.
denscomp(list(mmumax_weibull, mmumax_gamma , mmumax_lnorm, mmumax_norm), legendtext = plot.legend)
qqcomp(list(mmumax_weibull, mmumax_gamma , mmumax_lnorm, mmumax_norm), legendtext = plot.legend)
cdfcomp(list(mmumax_weibull, mmumax_gamma , mmumax_lnorm, mmumax_norm), legendtext = plot.legend)
ppcomp(list(mmumax_weibull, mmumax_gamma , mmumax_lnorm, mmumax_norm), legendtext = plot.legend)
par(mfrow = c(1, 1))

# Comparing different distributions - goodness-of-fit statistics.
gofstat(list(mmumax_weibull, mmumax_gamma , mmumax_lnorm, mmumax_norm), fitnames = c("weibull", "gamma", "lnorm", "norm"))
##-------------------------------------------------------
# sdMumax (actually it is ln(sdMumax))
plotdist(sdmumax_MCMCResults, histo = TRUE, demp = TRUE)
descdist(sdmumax_MCMCResults, boot = 1000)
# Fit candidate parametric distributions: Weibull, Gamma, Lognormal, and Normal.
sdmumax_weibull <- fitdist(sdmumax_MCMCResults, "weibull")
sdmumax_gamma <- fitdist(sdmumax_MCMCResults, "gamma")
sdmumax_lnorm <- fitdist(sdmumax_MCMCResults, "lnorm")
sdmumax_norm <- fitdist(sdmumax_MCMCResults, "norm")

# Comparing different distributions - goodness-of-fit plots.
par(mfrow = c(2, 2))
plot.legend <- c("Weibull", "gamma", "lognormal", "norm")

# Red---Weibull; Green---gamma; Blue---lognormal; Light blue---normal.
denscomp(list(sdmumax_weibull, sdmumax_gamma , sdmumax_lnorm, sdmumax_norm), legendtext = plot.legend)
qqcomp(list(sdmumax_weibull, sdmumax_gamma , sdmumax_lnorm, sdmumax_norm), legendtext = plot.legend)
cdfcomp(list(sdmumax_weibull, sdmumax_gamma , sdmumax_lnorm, sdmumax_norm), legendtext = plot.legend)
ppcomp(list(sdmumax_weibull, sdmumax_gamma , sdmumax_lnorm, sdmumax_norm), legendtext = plot.legend)
par(mfrow = c(1, 1))

# Comparing different distributions - goodness-of-fit statistics.
gofstat(list(sdmumax_weibull, sdmumax_gamma , sdmumax_lnorm, sdmumax_norm), fitnames = c("weibull", "gamma", "lnorm", "norm"))
##-------------------------------------------------------


