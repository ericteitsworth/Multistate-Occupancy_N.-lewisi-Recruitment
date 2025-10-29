#Teitsworth et al. 


#Set working directory
setwd("*")

#First, read in data
sitedata<-read.csv("SiteData.csv", header=TRUE)
obsdata<-read.csv("ObsData.csv", header=TRUE)
dethist<-read.csv("DetectionHistory.csv", header=TRUE)


# Following the Kery and Royle 2021 (ch6) code to get my data formatted and modeled 
# However, we aggregated all 5 years into "single season" multistate model to evaluate 
# hypotheses for the drivers of NRWD recruitment


##Summary of Data
#nsites = 176
#nyears = 5
#nsurveys = 4

#state 0 = unoccupied
#state 1 = occupied
#state 2 = occupied w/ evidence of reproduction


#Detection History
detdata<- dethist[,c(2:21)]

detdata <- detdata + 1 # Renumber observed states: 1,2,3 #Need to do this for JAGS
# Originally given as 0,1,2

# Put detection data and survey dates into a 3D array
(nsites <- nrow(detdata))
(nsurveys <- ncol(detdata))


library(AHMbook)

## Name covariates of interest
# Covs of psi, r
summary(ts<-sitedata$Total)
summary(sub<-sitedata$Bottom.Substrate)
summary(cover<-sitedata$Cover.Score)
summary(can<-sitedata$Light.Penetration) 
summary(pv<-sitedata$Pool.Variety)  
summary(cm<-sitedata$Channel.Modification..scaled.)
summary(bs<-sitedata$Bank.Stability...Vegetation..scaled.)
summary(rz<-sitedata$Riparian.Zone)
summary(tq<- dat2$TQ_AVG)

# Scale where needed: Standardize tool (AHMbook) both centers and scales, by default
ts.scaled<- standardize(ts)
sub.scaled<- standardize(sub)
cover.scaled<- standardize(cover)
can.scaled<- standardize(can)
tq.scaled<- standardize(tq)
pv.scaled<- standardize(pv)
rz.scaled<- standardize(rz)


#Detection covs
#p2, p32, p33
summary(discharge<-obsdata$Discharge.m)
discharge2<-discharge^2
summary(bait<-obsdata$Bait)
discharge.scaled<-standardize(discharge)
d2.scaled<-standardize(discharge2)
bait.scaled<-standardize(bait)


#Create an array of observation covariate data that is in the correct format
# to bundle in bdata below
d <- array(discharge.scaled, dim = c(nsites, nsurveys), dimnames = list(1:nsites, 1:nsurveys)) # this is nsites x nsurveys
d2<- array(d2.scaled, dim = c(nsites, nsurveys), dimnames = list(1:nsites, 1:nsurveys))
b <- array(bait.scaled, dim=c(nsites, nsurveys), dimnames = list(1:nsites, 1:nsurveys))


# Bundle data
str(bdata <- list(y = detdata, nsites = nsites, nsurveys = nsurveys,
                  sub=sub.scaled, cover=cover.scaled, can=can.scaled,rz=rz.scaled, bs=bs,
                  tq=tq.scaled,pv=pv.scaled, cm=cm,
                  d=d, d2=d2, b=b))



#### Model
# In this model, we are including the same covs in both psi and r

cat(file = "model.txt", "
model {

### (1) Linear Models and Priors
#-------------------------------
## (a) Latent State Linear Predictors

for (i in 1:nsites){
    logit(psi[i]) <- alpha.lpsi + beta.lpsi[1] * sub[i] + beta.lpsi[2] * cover[i] + beta.lpsi[3] * tq[i] + beta.lpsi[4] * cm[i] + beta.lpsi[5] * can[i] 
                    + beta.lpsi[6] * rz[i] + beta.lpsi[7] * pv[i] + beta.lpsi[8] * bs[i]
    logit(r[i]) <- alpha.lr + beta.lr[1] * sub[i] + beta.lr[2] * cover[i] + beta.lr[3] * tq[i] + beta.lr[4] * cm[i] + beta.lr[5] * can[i] 
                    + beta.lr[6] * rz[i] + beta.lr[7] * pv[i] + beta.lr[8] * bs[i]
}


# Priors for parameters in the linears models
    mean.psi ~ dunif(0, 1)           #occupancy intercept on probability scale
    alpha.lpsi <- logit(mean.psi)    #occupancy intercept
    mean.r ~ dunif(0, 1)             #repcruitment intercept on probability scale
    alpha.lr <- logit(mean.r)        #recruitment intercept
  
# Priors for covariates of parameters
  for(d in 1:8){
    beta.lpsi[d] ~ dnorm(0, 0.16)
  }  
  for(k in 1:8){  
    beta.lr[k] ~ dnorm(0, 0.16)
  }




## (b) Observation process
#Observation matrix Theta
#Linear models for p
for (i in 1:nsites){
  for (j in 1:nsurveys){
   #Observation model for sites in occupied state (state 2)
   logit(p2[i,j]) <- alpha.lp2 + beta.lp2[1] * d[i,j] + beta.lp2[2] * d2[i,j] + beta.lp2[3] * b[i,j]
    #Observation model for sites in recruiting state (state 3)
    mlogit.p3[2,i,j] <- alpha.lp32 + beta.lp32[1] * d[i,j] + beta.lp32[2] * d2[i,j] + beta.lp32[3] * b[i,j]          
    mlogit.p3[3,i,j] <- alpha.lp33 + beta.lp33[1] * d[i,j] + beta.lp33[2] * d2[i,j] + beta.lp33[3] * b[i,j]
  }
}                                                                                                                 
  
# Priors  
    mean.p2 ~ dunif(0, 1) #detection of true state 2 intercept on probability scale
    alpha.lp2 <- logit(mean.p2)   #detection of true state 2
    alpha.lp32 ~ dnorm(0, 0.1)    #Must be normal for multinomial logit
    alpha.lp33 ~ dnorm(0, 0.1)

  

  #Effect of discharge (continuous) and bait
  for(k in 1:3){                
    beta.lp2[k] ~ dnorm(0, 0.16)
    beta.lp32[k] ~ dnorm(0, 0.16)
    beta.lp33[k] ~ dnorm(0, 0.16)
  }
 
# Implement Multinomial logit link for p3[2:3] in Theta
  for (i in 1:nsites){
    for (j in 1:nsurveys){   # nsurveys is a function of site and year
      p3[2,i,j] <- exp(mlogit.p3[2,i,j]) / (1 + exp(mlogit.p3[2,i,j]) + exp(mlogit.p3[3,i,j]))
      p3[3,i,j] <- exp(mlogit.p3[3,i,j]) / (1 + exp(mlogit.p3[2,i,j]) + exp(mlogit.p3[3,i,j]))
    }
  }



### (2) Definte relationships between basic model structure and parameters
#-------------------------------------------------------------------------    
# Define state vector (Omega)
for(i in 1:nsites){
  Omega[i,1] <- 1 - psi[i]        # Prob. of non-occupation
  Omega[i,2] <- psi[i] * (1-r[i]) # Prob. of occupancy (w/ adults only)
  Omega[i,3] <- psi[i] * r[i]     # Prob. of occupancy (with recruiting)
}

# Define observation matrix (Theta)
# Order of indices: site, occasion, true state, observed state,
for(i in 1:nsites){
  for(j in 1:nsurveys){
    Theta[i,j,1,1] <- 1
    Theta[i,j,1,2] <- 0
    Theta[i,j,1,3] <- 0
    Theta[i,j,2,1] <- 1-p2[i,j]
    Theta[i,j,2,2] <- p2[i,j]
    Theta[i,j,2,3] <- 0
    Theta[i,j,3,1] <- 1-p3[2,i,j]-p3[3,i,j]
    Theta[i,j,3,2] <- p3[2,i,j]                                                            
    Theta[i,j,3,3] <- p3[3,i,j]
  }
}



### (3) Likelihoods
#-------------------------------------------------------
# State equation: model of true states (z)
for (i in 1:nsites){
  z[i] ~ dcat(Omega[i,])
}


# Observation
  for(i in 1:nsites){
    for(j in 1:nsurveys){
      y[i,j] ~ dcat(Theta[i,j,z[i],])
    }
  }
  


### (4) Derived quantities
# Number of sites in each state, conditional on observed data
for (i in 1:nsites){
  occ1[i] <- equals(z[i], 1)
  occ2[i] <- equals(z[i], 2)
  occ3[i] <- equals(z[i], 3)
  
}
  n.occ[1] <- sum(occ1[]) # Sites in state 1
  n.occ[2] <- sum(occ2[]) # Sites in state 2
  n.occ[3] <- sum(occ3[]) # Sites in state 3

# Detection probability of p32 and p33 on the probability scale
p3_derived_2 <- exp(alpha.lp32)/(1+exp(alpha.lp32) + exp(alpha.lp33))
p3_derived_3 <- exp(alpha.lp33)/(1+exp(alpha.lp32) + exp(alpha.lp33))

}

")

# Initial values
zst <- rep(3, nrow(bdata$y)) # Initialize at highest possible state
inits <- function(){list(z = zst)}

# Parameters monitored (could add "z")
params <- c("mean.psi","beta.lpsi[1]", "beta.lpsi[2]","beta.lpsi[3]","beta.lpsi[4]","beta.lpsi[5]","beta.lpsi[6]","beta.lpsi[7]","beta.lpsi[8]",
            "mean.r","beta.lr[1]", "beta.lr[2]","beta.lr[3]","beta.lr[4]","beta.lr[5]","beta.lr[6]","beta.lr[7]","beta.lr[8]",
            "mean.p2", "alpha.lp32", "alpha.lp33",
            "beta.lp2[1]", "beta.lp2[2]", "beta.lp2[3]",
            "beta.lp32[1]", "beta.lp32[2]", "beta.lp32[3]",
            "beta.lp33[1]", "beta.lp33[2]", "beta.lp33[3]",
            "n.occ", "occ1", "occ2", "occ3", "p3_derived_2","p3_derived_3"
)

# MCMC settings
na <- 1000 ; ni <- 150000 ; nt <- 20 ; nb <- 50000 ; nc <- 3  
# Call JAGS, check convergence and summarize posteriors # ~24hrs with 3 cores (intel 3.80GHz) and 32GB RAM 
library(jagsUI)
model <- jags(bdata, inits, params, "model.txt", n.adapt = na,
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

write.csv(model$summary, "*/model.csv")

