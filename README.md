# Uncovering latent territorial structure in ICFES Saber 11 performance with Bayesian multilevel spatial models

This repository contains a fully open source workflow using Markov chain Monte Carlo methods developed for *Uncovering latent territorial structure in ICFES Saber 11 performance with Bayesian multilevel spatial models*.

## Summary 

<p align="justify">
This article develops a Bayesian hierarchical framework to analyze academic performance in the 2022 second semester Saber 11 examination in Colombia. Our approach combines multilevel regression with municipal and departmental spatial random effects, and it incorporates Ridge and Lasso regularization priors to compare the contribution of sociodemographic covariates. Inference is implemented in a fully open source workflow using Markov chain Monte Carlo methods, and model behavior is assessed through synthetic data that mirror key features of the observed data. Simulation results indicate that Ridge provides the most balanced performance in parameter recovery, predictive accuracy, and sampling efficiency, while Lasso shows weaker fit and posterior stability, with gains in predictive accuracy under stronger multicollinearity. In the application, posterior rankings show a strong centralization of performance, with higher scores in central departments and lower scores in peripheral territories, and the strongest correlates of scores are student level living conditions, maternal education, access to educational resources, gender, and ethnic background, while spatial random effects capture residual regional disparities. A hybrid Bayesian segmentation based on K means propagates posterior uncertainty into clustering at departmental, municipal, and spatial scales, revealing multiscale territorial patterns consistent with structural inequalities and informing territorial targeting in education policy.
</p>

## Contents
<p align="justify">
This project contains the ICFES Saber 11 2022-II training and test datasets, three scenario datasets, and Colombia shapefiles at the departmental and municipal levels, which constitute the starting point of the analysis. The code used for data processing is also available; however, the files are already clean and prepared. Hybrid Metropolis–Hastings within Gibbs sampler codes developed for each model are provided, along with the code used to assess model fitting through simulations and real data applications. Although the code used to build the simulated datasets and spatial structures is available, all necessary objects and information are already included.

The following is the workflow to reproduce the results obtained in the analysis.
</p>

- Preparation
  - Load `data_ready.RData`,`MGN_ADM_MPIO_GRAFICO.shp`,`Spatial_Component.RData` and `MGN_DPTO_POLITICO.shp`.
  - Run `Descriptive_statistics_maps.R`.

- Samplers
  - For the ICFES Saber 11 training dataset and the three scenario datasets, run `G_sampler_M1_v6_C++.R` to obtain the MCMC output for the first model. The `Samplers.cpp` file must be uploaded beforehand. The dataset and the MCMC settings, such as the burn-in period and thinning step, are specified as arguments within the code.
  - For the Ridge model, run `G_sampler_Ridge_last_version.R`. The code uses the same arguments mentioned above and should be executed for the ICFES Saber 11 training dataset and the three scenario datasets.
  - For the Lasso model, run `G_sampler_Lasso_v2.R`. The code uses the same arguments mentioned above and should be executed for the ICFES Saber 11 training dataset and the three scenario datasets.
    
- Simulation assessment
  - Run `Convergence_diagnosis_M2_simulation.R` to obtain the baseline's evluation metrics on the three simulated dataset.
  - Run `Convergence_Ridge.R` to obtain the Ridge's evluation metrics on the three simulated dataset.
  - Run `Convergence_diagnosis_Laso_simulation.R` to obtain the Lasso's evluation metrics on the three simulated dataset.

- Real data assessment
  - For the goddness of fit, and predictive perfomance evaluation run `Convergence_Evaluation_all_models.R` and `Posterior_inference_analysis.R`.
  - For the inference analysis run `Results_Real_data_Rige.R`.


## Requirements
- This projec was developed on R / C++ (version R.4.5.1.).

