## ------------------------------------------------------------------------
#  Recall risk model (for LM in CSS)
## ------------------------------------------------------------------------
## Script purpose: Baseline model
## ------------------------------------------------------------------------
library(tidyverse);library(truncnorm);library(lhs);library(ggpubr);library(mc2d);library(Dict)
library(sensitivity);library(epiR);library(survival);library(Matching);library(quanteda);library(NCmisc)
## ------------------------------------------------------------------------
cat(strftime(Sys.time(),"%Y-%m-%d %H:%M:%S", tz = ""),'\n'," package import",'\n')
## ------------------------------------------------------------------------
# Define primary growth models and die-off & regrowth models

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
## Define other equations
# Function for converting Mumax from reference temperature: 25C
muAtNewTemp_CSS <- function(newTemp, oldMu, oldTemp = 25, T0 = -2.86) {
  numerator <- newTemp - T0
  denom <- oldTemp - T0
  newMu <- ((numerator / denom)^2) * oldMu
  return(newMu)
}
# Random sampling from truncated weibull distribution (for sampling LOG10Nmax for different nisin concentrations).
rtweibull <- function(n, shape, scale, a=0, b=Inf)
{
  stopifnot(n > 0 & all(scale > 0) & all(shape > 0))
  x <- runif(n)
  Fa <- pweibull(a, shape, scale)
  Fb <- pweibull(b, shape, scale)
  y <- (1-x)*Fa + x*Fb
  return(qweibull(y, shape, scale))
}

# Getting quantiles of truncated weibull distribution (for use in the sensitivity analysis).
qtweibull <- function(p, shape, scale, a=0, b=Inf)
{
  stopifnot( all(p >= 0 & p <= 1) & all( scale > 0 ) & all(shape > 0) )
  Fa <- pweibull(a, shape, scale)
  Fb <- pweibull(b, shape, scale)
  pNew <- p * (Fb - Fa) + Fa
  x <- qweibull(pNew, shape, scale)
  return(x)
}

# LHS sampling with restriction that within each sampling all probabilities add up to 1.
qdirichlet <- function(X, alpha) 
{ 
  # qdirichlet is not an exact quantile function since the quantile of a 
  #  multivariate distribtion is not unique 
  # qdirichlet is also not the quantiles of the marginal distributions since 
  #  those quantiles do not sum to one 
  # qdirichlet is the quantile of the underlying gamma functions, normalized 
  # This has been tested to show that qdirichlet approximates the dirichlet 
  #  distribution well and creates the correct marginal means and variances 
  #  when using a latin hypercube sample 
  lena <- length(alpha) 
  stopifnot(is.matrix(X)) 
  sims <- dim(X)[1] 
  stopifnot(dim(X)[2] == lena) 
  if (any(is.na(alpha)) || any(is.na(X))) 
    stop("NA values not allowed in qdirichlet") 
  
  Y <- matrix(0, nrow = sims, ncol = lena) 
  ind <- which(alpha != 0) 
  for (i in ind) 
  { 
    Y[,i] <- qgamma(X[,i], alpha[i], 1) 
  } 
  Y <- Y / rowSums(Y) 
  return(Y) 
} 
cat(strftime(Sys.time(),"%Y-%m-%d %H:%M:%S", tz = ""),'\n',"functions defined",'\n')
## ------------------------------------------------------------------------

### The baseline model ###

## ------------------------------------------------------------------------
## Set seed for reproducibility.
set.seed(05042020)
## ------------------------------------------------------------------------
## Prepare model-associated parameter values and distributions.

# Growth parameter dict.
growth_parameter <- read.delim("growth_parameter_0ppm")
cat(strftime(Sys.time(),"%Y-%m-%d %H:%M:%S", tz = ""),'\n',"parameter files read in",'\n')

# Strain dict.
serotype_vec <- c("12a", "12b", "4b")
strain_dict <- dict("0" = c("c1-111", "f2-237", "f6-366", "l3-051"), "5" = c("c1-111", "f6-366", "l3-051"),
                    "10" = c("c1-111", "f2-237", "f6-366", "l3-051"), "20" = c("c1-111", "f2-237", "f6-366", "l3-051"))

# Maximum population density.
LOG10Nmax <- dict("0" = dict(shape=41.341902, scale=9.126792), "5" = dict(shape=21.59984, scale=8.712891),
                  "10" = dict(shape=14.54644, scale=8.497603), "20" = dict(shape=8.438854, scale=8.208065))

# Product-associated parameters.
lot_iter_number=1000;package_number=10000;
net_wt=100;nisin_concentration=0;
fac_day=10;retail_day=30;house_day=20;

# Sampling-associated parameters.
samp_iter_number=10000;sampling_size=25;sampling_number=10;
sample_r=5;fac_freq=sample_r^2/(sample_r^2*fac_day+sample_r*retail_day+house_day);retail_freq=fac_freq/sample_r;house_freq=fac_freq/sample_r^2

# Prevalence.
prev <- 0.04016949; package_number_contam <- ceiling(package_number * prev)

# Temperature.
mtemp <- 4.4; sdtemp <- 0.77

# Minimum growth temperature.
mTmin <- -2.86; ln_sdTmin <- 0.638; sdTmin <- exp(ln_sdTmin)

# Initial contamination level.
mLnN0 <- -1.779525; sdLnN0 <- 1.068108

# Maximum growth rate.
mmumax <- 6.7091025; ln_sdmumax <- 0.9331415; sdmumax <- exp(ln_sdmumax)

cat(strftime(Sys.time(),"%Y-%m-%d %H:%M:%S", tz = ""),'\n',"parameters predefined",'\n')
## ------------------------------------------------------------------------
## Start of the model.
set.seed(05042020)

# Assign parameters to all packages in a lot-by-lot manner.
LOG10Nmax_list <- list()
temp_list <- list()
Tmin_list <- list()
N0_list <- list()
mumax_list <- list()
serotype_list <- list()
contam_list <- list()
sampling_day_series_list <- list()

for (lot in 1:lot_iter_number) {
  LOG10Nmax_list[[lot]] <- rtweibull(package_number, shape = LOG10Nmax$get(as.character(nisin_concentration))$get("shape"), scale = LOG10Nmax$get(as.character(nisin_concentration))$get("scale"), a=0, b=9.5)
}

for (lot in 1:lot_iter_number) {
  temp_list[[lot]] <- rnorm(package_number, mean = mtemp, sd = sdtemp)
  Tmin_list[[lot]] <- rnorm(package_number, mean = mTmin, sd = sdTmin)
  N0_list[[lot]] <- exp(rtruncnorm(package_number, a=log(1/net_wt), b=Inf, mean = mLnN0, sd = sdLnN0))
  mumax_list[[lot]] <- rtruncnorm(package_number, mean = mmumax, sd = sdmumax, a=0, b=Inf)/log(10)
  serotype_list[[lot]] <- sample(x=serotype_vec, package_number, replace = T)
  contam_list[[lot]] <- sample(x=c(1:package_number), package_number_contam, replace = F)
  sampling_day_series_list[[lot]] <- sample(x= c(1:(fac_day+retail_day+house_day)),size= samp_iter_number, 
                                            replace = TRUE, prob = c(rep(fac_freq,fac_day), rep(retail_freq, retail_day), rep(house_freq, house_day)))
}
cat(strftime(Sys.time(),"%Y-%m-%d %H:%M:%S", tz = ""),'\n',"parameters assigned to each lot",'\n') 

# Create vectors to store the results (expected values and credible intervals).
cumul_recall_risk <- rep(NA, lot_iter_number)
lower_CI <- rep(NA, lot_iter_number)
upper_CI <- rep(NA, lot_iter_number)

Baseline_Results_Dict <- dict(Cumulative_Recall_Risk = dict(Description = "Recall risk for each lot"), Lot_Test_Results = dict(Description = "Summary of test results for each lot"),
                  Package_Information = dict(Description = "Summary of package information for each lot"))

cat(strftime(Sys.time(),"%Y-%m-%d %H:%M:%S", tz = ""),'\n',"Output Object created",'\n')


# Start of the model.
seed = 05042020

for (lot in 1:lot_iter_number) {
  
  # Assign a unique seed for each lot.
  set.seed(seed+lot)
  
  # For each sampling, randomly pick a subset of samples to be tested.
  sample_ind_list <- list()
  for (sampling in 1:samp_iter_number) {
    sample_ind_list[[sampling]] <- sample(x=c(1:package_number), sampling_number, replace = F)
  }
  
  
  ## Inside lot loop.
  
  # Assign parameters for this lot (15 columns).
  package_df <- data.frame("Package_number" = c(1:package_number),
                           "Contamination" = rep(NA, package_number),
                           "Temperature" = temp_list[[lot]],
                           "Tmin" = Tmin_list[[lot]],
                           "N0" = N0_list[[lot]],
                           "mumax" = mumax_list[[lot]],
                           "new_mumax" = rep(NA, package_number),
                           "LOG10Nmax" = LOG10Nmax_list[[lot]],
                           "Serotype" = serotype_list[[lot]],
                           "Strain" = rep(NA, package_number), "growth_model" = rep(NA, package_number),
                           "p" = rep(NA, package_number), "delta" = rep(NA, package_number), 
                           "lag" = rep(NA, package_number),
                           "tc" = rep(NA, package_number))
  
  
  # Assign the contaminated packages to this lot.
  contam_package <- contam_list[[lot]]
  
  # In to package loop.
  ## ------------------------------------------------------------------------
  for (pack in 1:nrow(package_df)) {
    
    # For package pack, determine whether it is contaminated.
    
    if (pack %in% contam_package) {
      
      # If this package is contaminated, assign N0 (CFU/g) as the initial contamination level.
      package_df$Contamination[pack] <- 'positive'
      
      # Assign a representative serotype to this package.
      package_df$Serotype[pack] <- serotype <- sample(x=serotype_vec, 1, replace = T)
      
      if (package_df$Serotype[pack] == "4b") {
        package_df$Strain[pack] <- "f6-366"
      } else if (package_df$Serotype[pack] == "12b") {
        package_df$Strain[pack] <- "l3-051"
      } else if (length(strain_dict$get(as.character(nisin_concentration)))==4) {
        package_df$Strain[pack] <- sample(x=c("f2-237", "c1-111"), 1, replace = T)
      } else {
        package_df$Strain[pack] <- "c1-111"
      }
      
      # Extract the growth parameters for this package according to the representative strain.
      growth_model_strain <- growth_parameter[which(growth_parameter$Strains == package_df$Strain[pack]),]
      
      if (nrow(growth_model_strain)==1) {
        
        # If there is only one primary growth model, assign parameters.
        package_df$growth_model[pack] <- as.character(growth_model_strain$Model[1])
        package_df$p[pack] <- growth_model_strain$p[1]
        package_df$delta[pack] <- growth_model_strain$delta[1]
        package_df$lag[pack] <- growth_model_strain$lag[1]
        package_df$tc[pack] <- growth_model_strain$tc[1]
        
      } else {
        # If there are two primary growth models, determine which model to go with and assign parameters.
        a <- growth_model_strain$Model_wt[1]
        b <- growth_model_strain$Model_wt[2]
        growth_model_x <- runif(1,0,a+b)
        
        if (growth_model_x<=a) {
          # If the first model is selected
          package_df$growth_model[pack] <- as.character(growth_model_strain$Model[1])
          package_df$p[pack] <- growth_model_strain$p[1]
          package_df$delta[pack] <- growth_model_strain$delta[1]
          package_df$lag[pack] <- growth_model_strain$lag[1]
          package_df$tc[pack] <- growth_model_strain$tc[1]
        } else {
          # If the second model is selected
          package_df$growth_model[pack] <- as.character(growth_model_strain$Model[2])
          package_df$p[pack] <- growth_model_strain$p[2]
          package_df$delta[pack] <- growth_model_strain$delta[2]
          package_df$lag[pack] <- growth_model_strain$lag[2]
          package_df$tc[pack] <- growth_model_strain$tc[2]
        }}
      
    } else {
      # If package pack is not contaminated, assign 0 (CFU/g) as the initial contamination level.
      package_df$Contamination[pack] <- 'negative'
      package_df$N0[pack] <- 0
      package_df$Serotype[pack] <- NA
      
    }
    
    
    # Convert mumax to new_mumax for package pack.
    
    if (package_df$Temperature[pack] < package_df$Tmin[pack]) {
      package_df$new_mumax[pack] <- 0
    } else {
      package_df$new_mumax[pack] <- muAtNewTemp_CSS(newTemp=package_df$Temperature[pack], oldMu=package_df$mumax[pack], oldTemp = 25, T0 = package_df$Tmin[pack])
    }
    
    
  }
  ## ------------------------------------------------------------------------
  # Out of package loop.
  
  
  # Inside lot loop.
  
  # Create a matrix to store test results for this lot.
  lot_test_results <- data.frame(matrix(nrow = package_number, ncol = (fac_day+retail_day+house_day)))
  for (pack in 1:nrow(lot_test_results)) {
    rownames(lot_test_results)[pack] <- paste("package_", pack, sep = "")
  }
  for (day in 1:ncol(lot_test_results)) {
    colnames(lot_test_results)[day] <- paste("day_",day,sep = "")
  }
  
  
  # In to shelf life loop.
  ## ------------------------------------------------------------------------
  
  # For each day in shelf life, calculate the LM contamination level for each of the packages.
  for (day in 1:(fac_day+retail_day+house_day)) {
    
    N_vec <- rep(NA, nrow(package_df))
    
    # Into package loop
    ## ------------------------------------------------------------------------
    for (pack in 1:nrow(package_df)) {
      
      if (package_df$Contamination[pack] == 'positive') {
        # If package pack is contaminated, calculate the LOG10N on this day.
        LOG10N <- LOG10N_calculation(model_name = as.character(package_df$growth_model[pack]), t=day, 
                                     p=package_df$p[pack], delta=package_df$delta[pack],lag=package_df$lag[pack],
                                     mumax=package_df$new_mumax[pack],LOG10N0=log10(package_df$N0[pack]),
                                     LOG10Nmax=package_df$LOG10Nmax[pack],tc=package_df$tc[pack]) #+ random_error_list[[lot]][pack,day]
        
        # Convert LOG10N to # of cells per package.
        N_vec[pack] <- 10^LOG10N*net_wt
        
      } else {
        
        # If package pack is not contaminated, assign 0 on this day.
        N_vec[pack] <- 0
        
      }}
    
    ## ------------------------------------------------------------------------
    # Out of package loop.
    
    # Inside shelf life loop.
    
    package_df$day_x <- N_vec
    
    # For day day, change the column name into "day_day".
    names(package_df)[names(package_df) == "day_x"] <- paste("N_day_",day,sep = "")
    
  }
  
  ## ------------------------------------------------------------------------
  # Out of shelf life loop.
  
  
  # Inside lot loop.
  
  
  # In to package loop.
  ## ------------------------------------------------------------------------
  # If the # of cells per package goes below 1 on a specific day, assign 0 cell to the package for
  # that day and the following days.
  for (pack in 1:nrow(package_df)) {
    if (package_df$Contamination[pack] == "positive") {
      
      # In to shelf life loop.
      ## ------------------------------------------------------------------------
      for (day in 1:(fac_day+retail_day+house_day)) {
        if (package_df[pack,day+15] < 1) {
          # For this pack, if on a specific day the number of LM cells falls under 1 cell/package,
          # set 0 (CFU/package) to all days from this day to the end of shelf life.
          package_df[pack,(day+15):ncol(package_df)] <- rep(0, ncol(package_df)-day-14)
          break
        } else {next}}
      
      ## ------------------------------------------------------------------------
      # Out of shelf life loop.
      
    } else {next}}
  ## ------------------------------------------------------------------------
  # Out of package loop.
  
  
  # Inside lot loop.
  
  
  # Calculate test results.
  # In to shelf life loop.
  ## ------------------------------------------------------------------------
  for (day in 1:(fac_day+retail_day+house_day)) {
    test_results_day <- rep(NA, nrow(package_df))
    
    # In to package loop.
    ## ------------------------------------------------------------------------
    for (pack in 1:nrow(package_df)) {
      
      # For this pack, take the largest integer below LM_number, assign it to N.
      N <- floor(package_df[,day+15][pack])
      
      # For this pack, set the test result to 'No'.
      test_results_day[pack] <- "No"
      
      # Assign the probability of being tested as positive for this pack.
      if (N < 22) {
        positive_prob <- 1-phyper(0, N, net_wt-N, sampling_size)
      } else {
        positive_prob <- 0.999
      }
      
      # Determine whether this package should be tested as positive.
      positive_prob_x <- runif(1, 0, 1)
      
      if (positive_prob_x <= positive_prob) {
        test_results_day[pack] <- "Yes"
      }
      
    }
    
    ## ------------------------------------------------------------------------
    # Out of package loop
    
    
    lot_test_results[,day] <- test_results_day
    
  }
  
  ## ------------------------------------------------------------------------
  # Out of shelf life loop
  
  
  # Inside lot loop.
  
  
  # Retrieve sampling days.
  fac_freq=sample_r^2/(sample_r^2*fac_day+sample_r*retail_day+house_day)
  retail_freq=fac_freq/sample_r
  house_freq=fac_freq/sample_r^2
  sampling_day_series <- sampling_day_series_list[[lot]]
  
  
  # For each sampling, determine whether a recall event is resulted or not.
  # Calculate the number of recalls that should be issued.
  
  recall_number <- 0
  sample_id <- 0
  
  # In to sampling loop.
  ## ------------------------------------------------------------------------
  for (sample_day in sampling_day_series) {
    if (dim(table(lot_test_results[,sample_day])) == 2) {
      contam_ind <- which(lot_test_results[,sample_day]=="Yes")
      sample_id <- sample_id + 1
      sample_ind <- sample_ind_list[[sample_id]]
      contam_detect <- Reduce(intersect, list(contam_ind, sample_ind))
      if (length(contam_detect) != 0) {
        recall_number <- recall_number + 1
      }
    }}
  ## ------------------------------------------------------------------------
  # Out of sampling loop.
  
  recall_risk <- recall_number/samp_iter_number
  
  # Inside lot loop.
  
  cumul_recall_risk[lot] <- recall_risk
  
  Baseline_Results_Dict["Cumulative_Recall_Risk"][paste("Lot", lot, sep = "_")] = recall_risk
  Baseline_Results_Dict["Lot_Test_Results"][paste("Lot", lot, sep = "_")] = lot_test_results
  Baseline_Results_Dict["Package_Information"][paste("Lot", lot, sep = "_")] = package_df
  
  cat(strftime(Sys.time(),"%Y-%m-%d %H:%M:%S", tz = ""),'\n',paste("Lot number", lot, "is completed", sep = " "),'\n')
  
}


Baseline_Results_Dict$add("CRR" = cumul_recall_risk) 
# End of the model.

# Save the data.
save.image(file = "RecallRiskModel_Baseline_RC032521.RData")
## ------------------------------------------------------------------------

