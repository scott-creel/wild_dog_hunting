rm(list = ls())
setwd("C:/Users/g23b661/Desktop/DD_analysis/Ecosystems_years_combined_DD_covs")
#setwd("C:/Users/EcologyGuest1/OneDrive - Montana State University/Desktop/Ecosystems_years_combined_DD_covs")

library(ggplot2)
library(cowplot)
library(tidyverse)
library(viridis)
library(broom)
library(cowplot)
library(janitor)
library(diffdf)
library(glmmTMB)
library(ggsci)
library(DHARMa)
library(effects)
library(emmeans)
library(ggeffects)
library(scales)
library(amt)
library(chron)
library(readxl)

# read csv files created by the DD_add_covs scripts
# these have covariates values for rows from the ecosystem/year 
# identified in the filename  and NAs for the covariate values in all other rows

DD_gke_2022 <- read.csv("DD_dogs_with_all_covs_gke_2022.csv")
DD_gke_2023 <- read.csv("DD_dogs_with_all_covs_gke_2023.csv")
DD_lve_2022 <- read.csv("DD_dogs_with_all_covs_lve_2022.csv")
DD_lve_2023 <- read.csv("DD_dogs_with_all_covs_lve_2023.csv")

# deal with differences in columns created by prior code and joins

#identify mismatching columns
compare_df_cols(DD_gke_2022, DD_gke_2023, DD_lve_2022, DD_lve_2023)

#make the dataframes rowbind-able
DD_gke_2022 <- DD_gke_2022 %>% dplyr::select(!c(Distance_from_Lion_1, Distance_from_Lion_2, 
                                                Distance_from_Lion_3, DR.Alt, Mean_VeSBA, Outcome2))
DD_lve_2022 <- DD_lve_2022 %>% dplyr::select(!c(Distance_from_Lion_1, Distance_from_Lion_2, 
                                                Distance_from_Lion_3, DR.Alt, Mean_VeSBA, Outcome2))
DD_gke_2023 <- DD_gke_2023 %>% dplyr::select(!c(DR.Alt.x, VeSBA, VeSBA.smoothed,
                                                Outcome.for.second.half.of.HUNT, DR.Alt.x))
DD_lve_2023 <- DD_lve_2023 %>% dplyr::select(!c(DR.Alt.x, VeSBA, VeSBA.smoothed,
                                                Outcome.for.second.half.of.HUNT, DR.Alt.x))
DD_gke_2022 <- rename(DD_gke_2022, Outcome.of.HUNT = Outcome)
DD_lve_2022 <- rename(DD_lve_2022, Outcome.of.HUNT = Outcome)

# confirm row-bindable
compare_df_cols(DD_gke_2022, DD_gke_2023, DD_lve_2022, DD_lve_2023)
compare_df_cols_same(DD_gke_2022, DD_gke_2023, DD_lve_2022, DD_lve_2023)

#delete rows that are not pertinent to that ecosystem-year

DD_gke_2022 <- DD_gke_2022 %>% filter(Ecosystem == 'Kafue' & Year.y == 2022)
DD_gke_2023 <- DD_gke_2023 %>% filter(Ecosystem == 'Kafue' & Year.y == 2023)
DD_lve_2022 <- DD_lve_2022 %>% filter(Ecosystem == 'Luangwa' & Year.y == 2022)
DD_lve_2023 <- DD_lve_2023 %>% filter(Ecosystem == 'Luangwa' & Year.y == 2023)

DD_all <- rbind(DD_gke_2022, DD_gke_2023, DD_lve_2022, DD_lve_2023)
write.csv(DD_all, file = 'DD_all.csv')

# plot some random subsamples to confirm pattern is always the same
#DD_sub <- DD_all[sample(nrow(DD_all), 50000), ]

DD_all_narm <- DD_all %>% drop_na(VeDBA.smoothed, prot_extract, prey_biomass_extract,
                                  lion_extract, AdYear, Pups)
write.csv(DD_all_narm, file = 'DD_all_narm.csv')

DD_all_clean <- DD_all_narm %>% filter(prey_biomass_extract <150000)
write.csv(DD_all_clean, file = 'DD_all_clean.csv')

ar(DD_all_clean$VeDBA)

mypal <- c('darkgoldenrod3','cornsilk4','orange','red3','darkolivegreen',
           'cyan4', 'blue', 'violet', 'dodgerblue2', 'gray60', 'chartreuse',
           'chocolate3')

p <- ggplot(data = DD_all_clean, aes(x = prey_biomass_extract, y = VeDBA.smoothed, 
                                     group = 1, colour = as.factor(Pack))) +
  geom_point(alpha = 0.4) +
  guides(colour='none', alpha = 'none') +
  theme_classic() +
  scale_colour_manual(values = mypal) +
  xlab(bquote('Prey Biomass ('*'kg'~ km^-2*')')) +
  ylab("VeDBA (g)")+
  theme(axis.text = element_text (size=16))+
  theme(axis.title = element_text (size=20))

ggsave("vedba_all.tiff", device ="tiff")

mod1 <- glmmTMB(data = DD_all_clean, family = Gamma(link = "inverse"), 
                start = list(beta=c(0.1,0.01,0.01,0.01,0.01,0.01,0.01)),  REML = TRUE, 
                na.action = 'na.omit', formula = VeDBA.smoothed ~ scale(lag(VeDBA.smoothed, n=1)) +
                  scale(lion_extract) + scale(prot_extract) + scale(prey_biomass_extract) + 
                  scale(AdYear) + scale(Pups))

summary(mod1)

# coefficient backtrans, inverse is default for gamma (log or identity do not improve fit)

coeffs <- fixef(mod1)
bt_coeff <- 1/coeffs$cond  # for inverse link
#bt_coef <- exp(coeffs$cond) # for log link

#very similar backtransformed effects for the two links
#fit is better with inverse


# GOF with simulate
sims.gamma <- simulate(object=mod1, nsim = 100)
sim.vedba.mean <- colMeans(sims.gamma)
plot(density(sim.vedba.mean), xlim=c(0,0.5), main = 'Simulated vs Empirical Mean VeDBA',
     xlab="VeDBA")
abline(v=mean(DD_all$VeDBA.smoothed), col='orange', lwd=3)
abline(v=mean(sim.vedba.mean), col = 'blue', lwd=3, lty=2)

# very nice fit

# GOF with DHARMa 
DD_simres <- simulateResiduals(mod1)
plot(DD_simres)

#passes the dispersion test, scale-location alos looks OK
