# MBMARD
The Multivariate Bayesian Mixture AutoRegressive Decomposition

The MBMARD is the multivariate extension of the Bayesian Mixture Auto-Regressive Decomposition discussed in [Granados et al., 2022](https://www.sciencedirect.com/science/article/pii/S0167947321002437).
Analogously to BMARD is a stationary time series decomposition into latent processes such that each component is specified as a second-order autoregressive process AR(2). 

The AR(2) mixture components are now potentially shared among several observed series. In the application, modeling the EEG signals of two subjects allows first to decompose the EEG signals into a few main oscillations and group them based on the standard frequency bands used in neuroscience for human EEG. Then, the mixture approach allows us to identify which oscillations are shared with higher strength on the observed EEG nodes. Since the subjects performed a visual task, the occipital side of the brain showed the highest strengths. 

The details of the Multivariate BMARD method and the application can be found in [Granados et al., 2024](https://www.sciencedirect.com/science/article/pii/S2452306224000042)

## Examples

The method's performance can be tested by using the code in this repository; to that end, we share code to replicate the simulations found in the paper.

### Simulation_examples

The folder contains two R scripts with different examples of the MBMARD method.

Run the following code to simulate a mixture of AR(2) processes and fit the MBMARD model. A simple graphical summary of the spectral matrix estimation and the Partial Directed Coherence results is computed after fitting the model. 

MBMARD_Simulation_AR2mixture.R

This code replicates the same analysis for a mixture of AR(12), MA(4) and AR(1) processes. 

MBMARD_Simulation_ARMAmixture.R


### CPPcode
The main Rcpp code of the BMARD method, the primary function that runs the MCMC, is called MBMARD().

### Rcode_auxiliary
Auxiliary R code for MCMC post-processing is on the script MBMARDchainssummaryCurvecomputation.R, which computes the mean and median curves and the mean of the parameters determining the optimal number of components of the multichain run. The folder also contains the code to compute the partial directed coherence and the squared coherence of the model. The PDC is computed in the direction from the latent sources to the observed series.

