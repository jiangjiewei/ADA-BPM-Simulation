# ADA-BPM

This package provides the code and data for implementing our method and reproducing the numerical results in the article.

## Descriptions

The `ADA-BPM` folder contains files for reproducing the simulation studies:

- `ADA-BPM.R`, `Care.R`, `CD-trace.R`, `gCoda.R`, `SPIEC-EASI.R`: functions to implement precision matrix estimation for compositional data
- `Summary_roc.R`: functions to summarize the results of support recovery
- `Main_simulation.R`: functions to carry out the main simulations, reproducing Tables 1--2, A.4, and A.5
- `ROC_curve.R`: functions to plot ROC curves for compositional data, reproducing Figure 1

The `ada-bpm-real_data` folder contains files for reproducing the real data analysis:

- `basic_functions.R`, `proximal_gradient.R`, `tune_proximal_gradient.R`: functions to implement composition estimation using the method of Cao, Zhang and Li (2020)
- `ADA-BPM.R`, `Care.R`, `CD-trace.R`, `gCoda.R`, `SPIEC-EASI.R`: functions to implement precision matrix estimation for compositional data
- `Summary_network.R`: functions to summarize the identified networks
- `Main_real_data.R`: functions to carry out the real data analysis, reproducing Table 3 and Figures 2 
- `Gut Microbiome.RData`: gut microbiome data from Wu et al. (2011)
- `composition_hat.RData`, `networks.RData`: some intermediate results for fast computation

## Workflows

Before running the code, make sure that all the required R packages are successfully installed.

To carry out the main simulations for compositional data, open the `simulation` folder and run the R script:

- `Main_simulation.R`

Three graph structures (band, hub, and random) are considered for generating the basis precision matrix. Change the settings for dimension from 50 to 800 to obtain the desired results.
Please modify the names of output files when saving the results for different dimensions.

To apply the methods to gut microbiome data, open the `real_data` folder and run the R script:

- `Main_real_data.R`

Run the code for different methods to obtain the desired results. Some intermediate results are stored in `composition_hat.RData` and `networks.RData` to accelerate the computation.
