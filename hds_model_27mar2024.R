##new triage script

rm(list = ls())
library(sf)
library(rjags)
library(lubridate)
# if model isn't running, check for missing data values
# model and data object identifier

##############################################################
## CHANGE THESE VALUES TO MATCH THE INPUTS AND OUTPUT NAMES ##
##############################################################
runname<-paste0("GKE_buffalo_",Sys.Date(),"_",hour(Sys.time()),minute(Sys.time())) 
# species of interest
#spp=c("bushbuck","duiker","grysbok","medant","oribi","reedbuck","steenbok") 
spp=c("buffalo")



#basics, may need changing depending on needs
# cutoff distance in kms
B=0.5

transects<-st_read(dsn = "kaf_prey.gpkg", layer='observations')


lamCovs<-read.csv("geospatial_assets/gee_outputs/GKE_segments_quads_covariates.csv")
lamCovs<-lamCovs[,-(which(names(lamCovs)==".geo"))]


# the observations contain a new segment id - the new merged segments ( to ensure length > 1km) closest in space to the lat-long coordinates of the observation which is some cases puts the observer closer to a different segment than the one recorded by the observer - and also a different newly merged segment - so using these new segments to understand effort (number of times a site was surveyed) will be incomplete.  The code below addresses this.
original_segments<-st_read(dsn = "kaf_prey.gpkg", layer='gke_segments_original')
obs_old_segments<-as.numeric(unlist(lapply(strsplit(transects$segment,"_"),"[[",2)))

obs_new_segments<-original_segments$new_segmentid[match(paste0("kafue_",transects$segment),original_segments$segmentid)]# enpoints of transects were given an ID like a segment but aren't actually a segment.  identify all points (usually transect end points) that are not actually segments - these observations will need to get assigned an actual segment.
obs_old_segments[is.na(obs_new_segments)&obs_old_segments!=1]<-obs_old_segments[is.na(obs_new_segments)&obs_old_segments!=1]-1  # if these are end points they wont be 1, and we can subtract 1 to put them on an old segment that is an actual segment (not an end point).

obs_new_segments<- original_segments$new_segmentid[match(paste0("kafue_",transects$transect,"_",obs_old_segments),original_segments$segmentid)] # now match up every observation segment (as recorded by the observer) with the new segment ID for the sake of recording effort.

## make an effort table that better accounts for actual surveys
obs_table<-table(obs_new_segments,paste0(transects$studyseason,transects$studyyear))


overall_success<-rowSums(obs_table>0)[match(lamCovs$new_segmentid, row.names(obs_table))] # how many times in the entire study period of interest was a segment surveyed, reordered to match order of lamCovs

# restrict to buffer width
# need to think of  a solution for missing distances and group sizes.......
transects<-transects[which(transects$perp_distance<(B*1000)),]


# basic parameters for hds model
M=1000 # 1000 for impala and puku.
nind=sum(transects$species %in% spp)
nz=M-nind

#compile data
#### change here (name of data)
str(hds_data <- list(  # change name here for the current species
  B=B, # cutoff distance for distance sampling, in km
  M=M, # total size of augmented observations
  nind=nind, # number of observerd herds
  nz=nz, # number of augments
  groupsize=c(transects$group_size[transects$species %in% spp],rep(NA,nz))-1, # observered size of mixed group - 1 (for negbin)
  mixgroupsize=c(transects$mixed_group_size[transects$species %in% spp],rep(NA,nz))-1, # observered size of mixed group - 1 (for negbin)
  site= c(match(transects$herd_new_segmentid[transects$species %in% spp], lamCovs$new_segmentid),rep(NA,nz)), # interger site id
  d=c(transects$perp_distance[transects$species %in% spp],rep(NA,nz))/1000, # perpendicular distance in km between herd and transect line
  y=c(rep(1,nind),rep(0,nz)), # observed or not - all augments are 0
  z=c(rep(1,nind),rep(NA,nz)), # 'real' or not - all observations are 1, all augments are NA
  site_area=lamCovs$viewshed_area, # size of viewshed at each site
  nsites=length(lamCovs$new_segmentid), # number of sites
  treecover=lamCovs$treecover/100, # most common habitat scored at the site
  burnfreq=lamCovs$annburnfreq,
  dist_to_river=log(as.numeric(lamCovs$distriver)/1000), # distance to river in log km
  dist_to_stream=log(as.numeric(lamCovs$diststream)/1000),
  site_effort=(lamCovs$viewshed_area)*as.numeric(overall_success[match(lamCovs$new_segmentid,names(overall_success))]) # product of segment length and times surveyed
))

save(hds_data,file=paste0(runname,"_data.RData"))

# model
cat("
model{
  #priors
  alpha0 ~ dunif(-10,10) # intercept of log(sigma)in a half-normal detection function
  alpha1 ~ dunif(-10,10) # effect of groupsize-1 on log(sigma)
  
  beta0 ~ dunif(-10,10) # intercept for density (ind/km) - reference habitat is open woodland
  beta1 ~ dunif(-10,10) # lambda covariate for distance to river
  beta2 ~ dunif(-10,10) # lambda covariate for tree cover proportion
  beta3 ~ dunif(-10,10) # lambda covariate for annual burn probability
  beta4 ~ dunif(-10,10) # lambda covariate for distance to stream

  #Data Augmentation parameter
  psi <- sum(offlambda[]) / M # this is the only place absolute abundance plays a role
  
  #overdispersion for mixed group size (for detection)
  lambda.mixgroup ~ dgamma(0.1,0.1) #hyperprior for mean herd size
  rmix ~ dunif(0.1,100) #hyperprior for overdispersion
  pmixgroup<-rmix/(rmix+lambda.mixgroup)
  
  #overdispersion for species group size (for individal density)
  lambda.group ~ dgamma(0.1,0.1) #hyperprior for mean herd size
  r ~ dunif(0.1,100) #hyperprior for overdispersion
  pgroup<-r/(r+lambda.group)
  
  # individual(herd)-level detection model and site of occurence
  for(i in 1:M){
    groupsize[i] ~ dnegbin(pgroup, r) # for estimating individual density, not for detection 
    z[i] ~ dbern(psi) # does the ind exist
    d[i] ~ dunif(0,B) # distance (continuous) is uniform dist
    mixgroupsize[i] ~ dnegbin(pmixgroup, rmix) # for estimating detection 
    site[i] ~ dcat(site.probs[1:nsites]) # pop distribution among sites
    #zg[i] <- z[i]*(1 + groupsize[i])   # for estimate total numbers 
    log(sigma[i]) <- alpha0 + alpha1*pow(mixgroupsize[i],0.5)
    p[i] <- exp(-d[i]*d[i]/(2*sigma[i]*sigma[i])) # half norm det based on dist and parms
  }

  #observation model for each herd
  for(i in 1:M){
    y[i] ~ dbern(p[i]*z[i])
  }
  
  #models for herd density and site probability
  for (s in 1:nsites){ # for each site
    log(lambda[s]) <- beta0 + beta1*dist_to_river[s] + beta2*treecover[s] + beta3*burnfreq[s] + beta4*dist_to_stream[s] 
    offlambda[s]<-lambda[s]*site_effort[s] # density adjusted for size of the area and number of times it was surveyed 
    site.probs[s] <- offlambda[s]/sum(offlambda[]) # for dcat, the probability a random herd occurs at each site 
  }
  meanherddensity<-mean(lambda[1:nsites])
 
}",fill=TRUE, file="DA_allsurveys_p.mixgroupsize_lam.treeriverburn.txt")


library(rjags)
library(jagsUI)
ni <- 11000 ; nb <- 1000 ; nt <- 100 ; nc <- 10

# missing data?
# sometimes distance or group size is missing - it might be better to assign a random distance for that species on its specific segment if not part of the detection model

inits <- function(){list(alpha0=-1, alpha1=1, beta0=-1)} # for smallant, initiate beta0 at -5 any missing values need inits
params <- c("alpha0", "alpha1", "beta0","beta1","beta2","beta3","beta4","psi","lambda.group","r","lambda.mixgroup","rmix","meanherddensity")

# Call JAGS
#### change here (data type)

out<-jags(hds_data, inits, params,  "DA_allsurveys_p.mixgroupsize_lam.treeriverburn.txt", n.thin=nt,n.chains=nc, n.burnin=nb,n.iter=ni,parallel=T) # change name of species
print(out, 3)

save(out,file=paste0(runname,"_hds_out.RData"))



