## ------------------------------------------------------------------------
#  Recall risk model (for LM in CSS)
## ------------------------------------------------------------------------
## Script purpose: Sensitivity analysis.
##------------------------------------------------------
library(plyr);library(tidyverse);library(truncnorm);library(lhs);library(ggpubr);library(ggrepel)
library(sensitivity);library(epiR);library(survival);library(Matching);library(quanteda)
##------------------------------------------------------
## Load results of running the model under different scenarios.
Baseline_Results <- read.csv("./ModelRunningResults/CRR_Results/Baseline_CRR_RC032821.csv")
Var_Temp_Results <-read.csv("./ModelRunningResults/CRR_Results/var_temp_CRR_RC033021.csv")
Var_Tmin_Results <-read.csv("./ModelRunningResults/CRR_Results/var_Tmin_CRR_RC033021.csv")
Var_N0_Results <-read.csv("./ModelRunningResults/CRR_Results/var_N0_CRR_RC032821.csv")
Var_Mumax_Results <-read.csv("./ModelRunningResults/CRR_Results/var_mumax_CRR_RC033021.csv")
Var_LOG10Nmax_Results <-read.csv("./ModelRunningResults/CRR_Results/var_LOG10Nmax_CRR_RC033121.csv")
Unc_Temp_Results <- read.csv("./ModelRunningResults/CRR_Results/unc_temp_CRR_RC033021.csv")
Unc_Tmin_Results <- read.csv("./ModelRunningResults/CRR_Results/unc_Tmin_CRR_RC033021.csv")
Unc_N0_Results <- read.csv("./ModelRunningResults/CRR_Results/unc_N0_CRR_RC033021.csv")
Unc_Mumax_Results <- read.csv("./ModelRunningResults/CRR_Results/unc_mumax_CRR_RC033021.csv")
Unc_Prev_Results <- read.csv("./ModelRunningResults/CRR_Results/unc_prev_CRR_RC033021.csv")
Unc_Sampler_Results <- read.csv("./ModelRunningResults/CRR_Results/unc_sampler_CRR_RC033121.csv")
Unc_SeroFreq_Results <- read.csv("./ModelRunningResults/CRR_Results/unc_SeroFreq_CRR_RC033021.csv")
##-----------------------------------------------------------------------------
# Baseline statistics.
range(Baseline_Results$Baseline_CRR) # 0.3183 0.3489
median(Baseline_Results$Baseline_CRR)
quantile(Baseline_Results$Baseline_CRR, probs = c(0.025, 0.975))["2.5%"]
quantile(Baseline_Results$Baseline_CRR, probs = c(0.025, 0.975))["97.5%"]
# 0.3341 (0.3249975, 0.3432)
##-----------------------------------------------------------------------------
## Part I: Variable Parameters
##-----------------------------------------------------------------------------
## Calculate the SRCC between each variable parameter and predicted recall risk.

# Temp.
Var_Temp_Corr1 <- cor.test(x=Var_Temp_Results$temp, y=Var_Temp_Results$var_temp_CRR, method = 'spearman')
# Tmin.
Var_Tmin_Corr1 <- cor.test(x=Var_Tmin_Results$Tmin, y=Var_Tmin_Results$var_Tmin_CRR, method = 'spearman')
# N0.
Var_N0_Corr1 <- cor.test(x=Var_N0_Results$ln_N0, y=Var_N0_Results$var_N0_CRR, method = 'spearman')
# Mumaxref.
Var_Mumax_Corr1 <- cor.test(x=Var_Mumax_Results$mumax, y=Var_Mumax_Results$var_mumax_CRR, method = 'spearman')
# LOG10Nmax.
Var_LOG10Nmax_Corr1 <- cor.test(x=Var_LOG10Nmax_Results$LOG10Nmax, y=Var_LOG10Nmax_Results$Var_LOG10Nmax_CRR, method = 'spearman')
##-----------------------------------------------------------------------------
## Part II: Uncertain Parameters
##-----------------------------------------------------------------------------
## Calculate the SRCC/PRCC between each uncertain parameter and predicted recall risk.

# mtemp and sdtemp
Unc_Temp_Results_Selected <- subset(Unc_Temp_Results, select = -c(Lot))
LHS_Unc_Temp_PRCC <- epi.prcc(Unc_Temp_Results_Selected,sided.test = 2)
# mTmin and sdTmin
Unc_Tmin_Results_Selected <- subset(Unc_Tmin_Results, select = -c(Lot))
LHS_Unc_Tmin_PRCC <- epi.prcc(Unc_Tmin_Results_Selected,sided.test = 2)
# mLnN0 and sdLnN0
Unc_N0_Results_Selected <- subset(Unc_N0_Results, select = -c(Lot))
LHS_Unc_N0_PRCC <- epi.prcc(Unc_N0_Results_Selected,sided.test = 2)
# mMumaxref and sdMumaxref
Unc_Mumax_Results_Selected <- subset(Unc_Mumax_Results, select = -c(Lot, ln_sdmumax))
LHS_Unc_Mumax_PRCC <- epi.prcc(Unc_Mumax_Results_Selected,sided.test = 2)
# Prev
Unc_Prev_Corr1 <- cor.test(x=Unc_Prev_Results$prev, y=Unc_Prev_Results$unc_prev_CRR, method = 'spearman')
# Sampling_R.
Unc_Sampler_Corr1 <- cor.test(x=Unc_Sampler_Results$Sample_r, y=Unc_Sampler_Results$Unc_Sampler_CRR, method = 'spearman')
# Prop_12a, Prop_12b, and Prop_4b
Unc_SeroFreq_Results_Selected <- subset(Unc_SeroFreq_Results, select = -c(Lot))
LHS_Unc_SeroFreq_PRCC <- epi.prcc(Unc_SeroFreq_Results_Selected,sided.test = 2)

##-----------------------------------------------------------------------------
## Part III: Calculate BH-corrected P-values for each parameter.
##-----------------------------------------------------------------------------
InputFile <- read.csv("Parameters_CorrelationCoefficients_RC040121.csv", header = TRUE)
InputFile$BH_P <- p.adjust(InputFile$Naive_P,method="BH")
##-----------------------------------------------------------------------------

