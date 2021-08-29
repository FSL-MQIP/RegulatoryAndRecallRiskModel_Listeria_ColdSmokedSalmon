## ------------------------------------------------------------------------
#  Recall risk model (for LM in CSS)
## ------------------------------------------------------------------------
## Script purpose: Model validation - main model.
## ------------------------------------------------------------------------

## Start of the model.
## ------------------------------------------------------------------------
set.seed(03042021)

# Get parameters for all packages in a lot-by-lot manner.
## ------------------------------------------------------------------------
LOG10Nmax_list <- list()
for (lot in 1:lot_iter_number) {
  LOG10Nmax_list[[lot]] <- rtweibull(package_number, shape = LOG10Nmax$get(as.character(nisin_concentration))$get("shape"), scale = LOG10Nmax$get(as.character(nisin_concentration))$get("scale"), a=0, b=9.5)
}

sigma_vec <- exp(rnorm(lot_iter_number, mean = mlnsigma, sd = sdlnsigma))
temp_list <- list()
Tmin_list <- list()
LOG10N0_list <- list()
mumax_list <- list()
serotype_list <- list()
random_error_list <- list()

for (lot in 1:lot_iter_number) {
  temp_list[[lot]] <- runif(package_number, min = temp - 0.5, max = temp + 0.5)
  Tmin_list[[lot]] <- rnorm(package_number, mean = mTmin[lot], sd = sdTmin[lot])
  #Tmin_list[[lot]] <- rep(mTmin, package_number)
  #LOG10N0_list[[lot]] <- rweibull(package_number, shape = shape_LOG10N0, scale = scale_LOG10N0)
  LOG10N0_list[[lot]] <- rnorm(package_number, mean = mLOG10N0, sd = sdLOG10N0)
  mumax_list[[lot]] <- rtruncnorm(package_number, mean = mmumax[lot], sd = sdmumax[lot], a=0, b=Inf)/log(10)
  serotype_list[[lot]] <- sample(x=serotype_vec, package_number, replace = T)
  random_error_list[[lot]] <- rnorm(package_number, mean = 0, sd = sigma_vec[lot])
}
## ------------------------------------------------------------------------

# Start simulation.
## ------------------------------------------------------------------------

simulated_df <- data.frame()
simulated_data <- vector()

for (lot in 1:lot_iter_number) {
  
  # Inside lot loop.
  
  package_df <- data.frame("Package_number" = c((package_number*(lot-1)+1):(package_number*(lot-1)+package_number)),
                           "Temperature" = temp_list[[lot]],
                           "Tmin" = Tmin_list[[lot]],
                           "LOG10N0" = LOG10N0_list[[lot]],
                           "mumax" = mumax_list[[lot]],
                           "new_mumax" = rep(NA, package_number),
                           "LOG10Nmax" = LOG10Nmax_list[[lot]],
                           "Serotype" = serotype_list[[lot]],
                           "Strain" = rep(NA, package_number), "growth_model" = rep(NA, package_number),
                           "p" = rep(NA, package_number), "delta" = rep(NA, package_number), 
                           "lag" = rep(NA, package_number),
                           "tc" = rep(NA, package_number), "LOG10N" = rep(NA, package_number),
                           "RandomError" = random_error_list[[lot]])
  # 15 columns.
  
  # Into package loop.
  ## ------------------------------------------------------------------------
  for (pack in 1:nrow(package_df)) {
    
    if (package_df$Serotype[pack] == "4b") {
      package_df$Strain[pack] <- "f6-366"
    } else if (package_df$Serotype[pack] == "12b") {
      package_df$Strain[pack] <- "l3-051"
    } else if (length(strain_dict$get(as.character(nisin_concentration)))==4) {
      package_df$Strain[pack] <- sample(x=c("f2-237", "c1-111"), 1, replace = T)
    } else {
      package_df$Strain[pack] <- "c1-111"
    }
    
    # For this package, extract growth parameters according to the representative strain.
    growth_model_strain <- growth_parameter[which(growth_parameter$Strains == package_df$Strain[pack]),]
    
    if (nrow(growth_model_strain)==1) {
      
      # If there is only one primary growth model
      package_df$growth_model[pack] <- as.character(growth_model_strain$Model[1])
      package_df$p[pack] <- growth_model_strain$p[1]
      package_df$delta[pack] <- growth_model_strain$delta[1]
      package_df$lag[pack] <- growth_model_strain$lag[1]
      package_df$tc[pack] <- growth_model_strain$tc[1]
      
    } else {
      # If there are two primary growth models
      # Determine, for package pack, which model to go with
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
    
    
    if (package_df$Temperature[pack] < package_df$Tmin[pack]) {
      package_df$new_mumax[pack] <- 0
    } else {
      package_df$new_mumax[pack] <- muAtNewTemp_CSS(newTemp=package_df$Temperature[pack], oldMu=package_df$mumax[pack], oldTemp = 25, T0 = package_df$Tmin[pack])
    }
    
    
  }
  ## ------------------------------------------------------------------------
  # Out of package loop.
  
  
  # In to package loop.
  ## ------------------------------------------------------------------------
  for (pack in 1:nrow(package_df)) {
    
    package_df$LOG10N[pack] <- LOG10N_calculation(model_name = as.character(package_df$growth_model[pack]), t=day, 
                                                  p=package_df$p[pack], delta=package_df$delta[pack],lag=package_df$lag[pack],
                                                  mumax=package_df$new_mumax[pack],LOG10N0=package_df$LOG10N0[pack],
                                                  LOG10Nmax=package_df$LOG10Nmax[pack],tc=package_df$tc[pack]) + package_df$RandomError[pack]
    
  }
  ## ------------------------------------------------------------------------
  # Out of package loop.
  
  
  # Inside lot loop.
  
  simulated_df <- rbind(simulated_df, package_df)
  
  # Inside validation loop.
  
  simulated_data <- c(simulated_data, package_df$LOG10N)
  
  
}

## ------------------------------------------------------------------------
# End of simulation.

## ------------------------------------------------------------------------
# End of the model.