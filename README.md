# wild_dog_hunting
R, Java and JAGS code for analysis of effects on hunting costs and benefits

Code for analysis of the effects of prey density, lion density, protection level, pack size and the number of dependent pups on energy expendiure, movement, hunting success and prey size, including
1.  Java script (used within Google Earth Engine) to extract environmental variables used as predictors in hierarchical distance sampling models to estimate the density of each prey species.
2.  R script used to merge raw transect data woth covariates and fit Bayesian hierarchical distance models of prey density.
3.  JAGS code for the hierarchical distance sampling models fit witin R using jagsUI, saved as ASCII (.txt).
4.  R scripts to compile data, fit and validate gamma regression models to test for effects on hunting costs and benefits (example using vectorial dynamic body acceleration as the dependent variable).
