# RecallRiskModel_Listeria_ColdSmokedSalmon

## Purpose
This repository data and codes associated with developing and studying the Recall Risk Model.

## File description
### RRM_GH_BaselineModel_RC082821.R
This is the script of Recall Risk Model with parameters settings specific to the baseline scenario. Running the models under alternative scenarios (e.g., those for scenario analysis) can be achieved by modifying this script.

### RRM_GH_SensitivityAnalysis_RC082821.R
This is the script used for performing sensitivity analysis of the Recall Risk Model for identification of variable and uncertain parameters that are significantly correlated with the recall risks.

### RRM_GH_ModelValidation_RC082821.R
This is the script used for performing model validation of the Recall Risk Model.

### RRM_GH_ModelValidation_MainModel_RC082821.R
This is the script of the die-off/growth dynamics part of the Recall Risk Model modified for use in model validation. This script is called in "RRM_GH_ModelValidation_RC082821.R".

### RRM_GH_GCFitting_RC082821.R
This is a template script used for fitting a given curve with different primary growth or die-off/regrowth models and comparing different models to determine the best-fit model.

### RRM_GH_DitributionFitting_RC082821.R
This is a template script used for fitting parameter data with different parametric distributions and comparing different distributions to determine the best-fit distribution.

### RRM_GH_ComputePosterior_RC082821.R
This is a template script used for updating prior distributions of a given parameter with additional data to posterior distributions following a Bayesian approach.

### RecallRiskModeling_Results
This is a folder containing the recall risks results of running the Recall Risk Model under different scenarios, including the baseline scenarios, the scenarios designed for sensitivity analysis, and the scenarios designed for scenario analysis.

### GrowthParameters
This is a folder containing files showing the die-off/growth parameter values of different representative strains and nisin concentrations. Files in this folder are called in the scripts for running the Recall Risk Model or doing model validation.

### ValidationData
This folder contains a file showing the experimental data used for validating the Recall Risk Model (Chen et al., 2020).

### GrowthCurveData
This folder contains files showing the growth curve data (Kang et al., 2014) used for determining the die-off/growth parameter values.

### NewMumaxData
This folder contains new data of the maximum growth rate of L. monocytogenes on cold-smoked salmon from three different studies . These data were used to update the distributions of mMumaxref and sdMumaxref (hyperparameters of Mumaxref).

### InitialContamLevelData
This folder contains a file showing the data of initial contamination level of L. monocytogenes on cold-smoked salmon collected by Afssa.
