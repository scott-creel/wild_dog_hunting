# wild_dog_hunting
R, Java and JAGS code for analysis of effects on hunting costs and benefits

Code for analysis of the effects of prey density, lion density, protection level, pack size and the number of dependent pups on energy expendiure, movement, hunting success and prey size, including
1.  Java script (used within Google Earth Engine) to extract environmental variables used as predictors in hierarchical distance sampling models to estimate the density of each prey species.
2.  R script used to merge prey species detections in raw transect data with covariates and fit Bayesian hierarchical distance sampling models of prey density.
3.  JAGS code for the hierarchical distance sampling models fit within R using jagsUI, saved as ASCII (.txt).
4.  R scripts to compile data, fit and assess GOF for (gamma and poisson) regression models to test for effects on hunting costs and benefits (example using vectorial dynamic body acceleration as the dependent variable).
