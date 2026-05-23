#Benthic Cover, Riparian Zones, and Flashy Urban Streams Limit Juvenile Recruitment of a Long-lived Aquatic Salamander
#Teitsworth et al. 

##---------------------------------------------------------
# Plotting Multistate Model Result
##---------------------------------------------------------

# Plot effects of model to better visualize in a caterpillar plot
library(mcmcplots)
# Inner quantile is 68% CI and Outer quantile is 95% CI

# r effects
par(mfcol = c(1,1), mar=c(4,9,2,2), cex.axis=2.0)
caterplot(model$samples,parms = c("beta.lr[1]", "beta.lr[2]","beta.lr[3]","beta.lr[4]","beta.lr[5]","beta.lr[6]",
                                  "beta.lr[7]","beta.lr[8]"), regex=NULL, random=NULL, reorder=FALSE, col = "red", pch=19, cex=2.0, collapse=TRUE, cat.shift=-0.1,
          labels = c("","","","","","","",""),lwd=c(1.0,3.0), style="plain")

# psi effects
caterplot(model$samples, parms = c("beta.lpsi[1]", "beta.lpsi[2]","beta.lpsi[3]","beta.lpsi[4]","beta.lpsi[5]","beta.lpsi[6]",
                                   "beta.lpsi[7]","beta.lpsi[8]"), add=TRUE, collapse=TRUE, reorder=FALSE, pch=19, cex=2.0, labels = c(
                                     "Substrate","Cover","TQavg","Channel","Canopy","Riparian","Pools","Bank"), col="black", lwd=c(1.0,3.0))
abline(v=0, col="black")

# Detection effects
DetPlots <- par(mfcol = c(3,1), mar=c(3,10,2,2), cex.main=1.5)
caterplot(model$samples, parms = c("beta.lp2[1]", "beta.lp2[2]", "beta.lp2[3]"), regex=NULL, random=NULL, collapse=TRUE, lwd=c(1.0,3.0), pch=19, cex=1.5,
          reorder=FALSE, labels = c("Discharge", "Discharge^2", "Bait Age"),style="plain", col="black", xaxt="n", val.lim = c(-6,6), cex.main=1.5)
title(main="Detecting Adults at Non-Recruiting")
abline(v=0, col="black")
caterplot(model$samples, parms = c("beta.lp32[1]", "beta.lp32[2]", "beta.lp32[3]"), regex=NULL, random=NULL, collapse=TRUE,lwd=c(1.0,3.0), pch=19, cex=1.5,
          reorder=FALSE, labels = c("Discharge", "Discharge^2", "Bait Age"), col="black", style="plain",xaxt="n",val.lim = c(-6,6), cex.main=1.5)
abline(v=0, col="black")
title(main="Detecting Adults at Recruiting")
caterplot(model$samples, parms = c("beta.lp33[1]", "beta.lp33[2]", "beta.lp33[3]"), regex=NULL, random=NULL, collapse=TRUE,lwd=c(1.0,3.0), pch=19, cex=1.5,
          reorder=FALSE, labels = c("Discharge", "Discharge^2", "Bait Age"), style="plain", col="black",val.lim = c(-6,6), cex.axis=1.5, cex.main=1.5)
abline(v=0, col="black")
title(main="Detecting Evidence of Recruitment (Juveniles)")


## --------------------------------------------------
## Prediction plots - expected probabilities as a function of covariates
## --------------------------------------------------  
  
  #grab posterior means
  tmp<- model$samples
  
  ## mcmc output as a matrix
  outmat <- as.matrix(model$samples)
  
  ## number of prediction points
  npred<-500
  
  pred.psi.sub<-pred.psi.cover<-pred.psi.pv<-pred.psi.bs<-pred.psi.can<-
    pred.psi.cm<-pred.psi.tq<-pred.psi.rz<-
    pred.r.cover<-pred.r.tq<-pred.r.pv<-pred.r.rz<-
    pred.r.can<-pred.r.sub<-pred.r.cm<-pred.r.bs<-matrix(NA, nrow(outmat), npred)
  
  # prediction values for covariates  
  sub0 <- seq(min(sub), max(sub), length.out=npred)
  cover0 <- seq(min(cover), max(cover), length.out=npred)
  cm0 <- seq(min(cm), max(cm), length.out=npred)
  pv0 <- seq(min(pv), max(pv), length.out=npred)
  rz0 <- seq(min(rz), max(rz), length.out=npred)
  tq0 <- seq(min(tq), max(tq), length.out=npred)
  bs0 <- seq(min(bs), max(bs), length.out=npred)
  can0 <- seq(min(can), max(can), length.out=npred)
  
  ## standardize covariates relative to those used in analyses
  subP <- standardize2match(sub0, sub)
  coverP <- standardize2match(cover0, cover)
  cmP <- standardize2match(cm0, cm)
  pvP <- standardize2match(pv0, pv)
  rzP <- standardize2match(rz0, rz)
  tqP <- standardize2match(tq0, tq)
  bsP <- standardize2match(bs0, bs)
  canP <- standardize2match(can0, can)
  
  ## Predict at the means of the other covariates (i.e., 0)
  ## Here, predict the probabilities at each MCMC iteration
  ## then summarize the posterior distribution of each prediction
  ## this can take a couple minutes.
  for(i in 1:nrow(outmat)){
    pred.psi.sub[i,] <- plogis(qlogis(outmat[i,"mean.psi"]) + outmat[i,"beta.lpsi[1]"] * subP)
    pred.psi.cover[i,] <- plogis(qlogis(outmat[i,"mean.psi"]) + outmat[i,"beta.lpsi[2]"] * coverP)
    pred.psi.cm[i,] <- plogis(qlogis(outmat[i,"mean.psi"]) + outmat[i,"beta.lpsi[4]"] * cmP)
    pred.psi.pv[i,] <- plogis(qlogis(outmat[i,"mean.psi"]) + outmat[i,"beta.lpsi[7]"] * pvP)
    pred.psi.rz[i,] <- plogis(qlogis(outmat[i,"mean.psi"]) + outmat[i,"beta.lpsi[6]"] * rzP)
    pred.psi.tq[i,] <- plogis(qlogis(outmat[i,"mean.psi"]) + outmat[i,"beta.lpsi[3]"] * tqP)
    pred.psi.bs[i,] <- plogis(qlogis(outmat[i,"mean.psi"]) + outmat[i,"beta.lpsi[8]"] * bsP)
    pred.psi.can[i,] <- plogis(qlogis(outmat[i,"mean.psi"]) + outmat[i,"beta.lpsi[5]"] * canP)
    pred.r.sub[i,] <- plogis(qlogis(outmat[i,"mean.r"]) + outmat[i,"beta.lr[1]"] * subP)
    pred.r.cover[i,] <- plogis(qlogis(outmat[i,"mean.r"]) + outmat[i,"beta.lr[2]"] * coverP)
    pred.r.cm[i,] <- plogis(qlogis(outmat[i,"mean.r"]) + outmat[i,"beta.lr[4]"] * cmP)
    pred.r.pv[i,] <- plogis(qlogis(outmat[i,"mean.r"]) + outmat[i,"beta.lr[7]"] * pvP)
    pred.r.rz[i,] <- plogis(qlogis(outmat[i,"mean.r"]) + outmat[i,"beta.lr[6]"] * rzP)
    pred.r.tq[i,] <- plogis(qlogis(outmat[i,"mean.r"]) + outmat[i,"beta.lr[3]"] * tqP)
    pred.r.bs[i,] <- plogis(qlogis(outmat[i,"mean.r"]) + outmat[i,"beta.lr[8]"] * bsP)
    pred.r.can[i,] <- plogis(qlogis(outmat[i,"mean.r"]) + outmat[i,"beta.lr[5]"] * canP) 
  }
  
  
  ## --------------------------------------------------
  #Plot psi effects
  
  par(mfrow = c(2,4), mar=c(5,5,5,2), cex.lab = 1.75, cex.axis = 1.5)
  matplot(sub0, apply(pred.psi.sub,2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0, 15), col="red",
          xlab = "Substrate Quality", ylab = "Occupancy", 
          #main = "Substrate Quality effect on Occupancy", 
          frame = FALSE)
  lines(seq(min(sub), max(sub), length.out=npred), apply(pred.psi.sub,2,function(x) quantile(x, 0.975)), col="red", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(sub), max(sub), length.out=npred), apply(pred.psi.sub,2,function(x) quantile(x, 0.025)), col="red", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(cover0, apply(pred.psi.cover,2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(5, 20), col="red",
          xlab = "Cover Quality", ylab = "",
          #main = "Cover Quality effect on Occupancy", 
          frame = FALSE)
  lines(seq(min(cover), max(cover), length.out=npred), apply(pred.psi.cover,2,function(x) quantile(x, 0.975)), col="red", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(cover), max(cover), length.out=npred), apply(pred.psi.cover,2,function(x) quantile(x, 0.025)), col="red", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(pv0, apply(pred.psi.pv,2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0, 10), col="red",
          xlab = "Pool Variety", ylab = "",
          # main = "Pool Variety effect on Occupancy", 
          frame = FALSE)
  lines(seq(min(pv), max(pv), length.out=npred), apply(pred.psi.pv,2,function(x) quantile(x, 0.975)), col="red", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(pv), max(pv), length.out=npred), apply(pred.psi.pv,2,function(x) quantile(x, 0.025)), col="red", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(bs0, apply(pred.psi.bs,2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(-8, 6),
          xlab = "Bank Stability (scaled)", ylab = "",
          # main = "Bank Stability effect on Occupancy", 
          frame = FALSE)
  lines(seq(min(bs), max(bs), length.out=npred), apply(pred.psi.bs,2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(bs), max(bs), length.out=npred), apply(pred.psi.bs,2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(can0,apply(pred.psi.can,2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0, 10),
          xlab = "Canopy Quality", ylab = "Occupancy",
          # main = "Canopy Quality effect on Occupancy", 
          frame = FALSE)
  lines(seq(min(can), max(can), length.out=npred), apply(pred.psi.can,2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(can), max(can), length.out=npred), apply(pred.psi.can,2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(cm0, apply(pred.psi.cm,2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(-15, 4.25), col="red",
          xlab = "Channel Modification (scaled)", ylab = "",
          # main = "Channel Modification effect on Occupancy", 
          frame = FALSE)
  lines(seq(min(cm), max(cm), length.out=npred), apply(pred.psi.cm,2,function(x) quantile(x, 0.975)), col="red", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(cm), max(cm), length.out=npred), apply(pred.psi.cm,2,function(x) quantile(x, 0.025)), col="red", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(tq0, apply(pred.psi.tq,2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0.15, 0.4),
          xlab = "TQavg (% Year Above Mean)", ylab = "",
          #main = "Stream Flashiness effect on Occupancy", 
          frame = FALSE)
  lines(seq(min(tq), max(tq), length.out=npred), apply(pred.psi.tq,2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(tq), max(tq), length.out=npred), apply(pred.psi.tq,2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(rz0, apply(pred.psi.rz,2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(5, 10), col="red",
          xlab = "Riparrian zone", ylab = "",
          frame = FALSE)
  lines(seq(min(rz), max(rz), length.out=npred), apply(pred.psi.rz,2,function(x) quantile(x, 0.975)), col="red", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(rz), max(rz), length.out=npred), apply(pred.psi.rz,2,function(x) quantile(x, 0.025)), col="red", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  ## --------------------------------------------------
  # Plot juvenile recruitment effects
  # IMPORTANT - this is the probability of juvenile recruitment
  # It is NOT the probability that adults and juveniles are present. 
  # The latter is psi*r
  
  par(mfrow = c(2,4), mar=c(5,5,5,2), cex.lab = 1.75, cex.axis = 1.5)
  matplot(sub0, apply(pred.r.sub,2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0, 15), col="red",
          xlab = "Substrate Quality", ylab = "Recruitment", 
          frame = FALSE)
  lines(seq(min(sub), max(sub), length.out=npred), apply(pred.r.sub,2,function(x) quantile(x, 0.975)), col="red", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(sub), max(sub), length.out=npred), apply(pred.r.sub,2,function(x) quantile(x, 0.025)), col="red", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(cover0, apply(pred.r.cover,2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(5, 20), col="black",
          xlab = "Cover Quality", ylab = "",
          frame = FALSE)
  lines(seq(min(cover), max(cover), length.out=npred), apply(pred.r.cover,2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(cover), max(cover), length.out=npred), apply(pred.r.cover,2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(pv0, apply(pred.r.pv,2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0, 10), col="black",
          xlab = "Pool Variety", ylab = "",
          frame = FALSE)
  lines(seq(min(pv), max(pv), length.out=npred), apply(pred.r.pv,2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(pv), max(pv), length.out=npred), apply(pred.r.pv,2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(bs0, apply(pred.r.bs,2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(-8, 6), col="red",
          xlab = "Bank Stability (scaled)", ylab = "",
          frame = FALSE)
  lines(seq(min(bs), max(bs), length.out=npred), apply(pred.r.bs,2,function(x) quantile(x, 0.975)), col="red", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(bs), max(bs), length.out=npred), apply(pred.r.bs,2,function(x) quantile(x, 0.025)), col="red", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(can0,apply(pred.r.can,2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0, 10), col="red",
          xlab = "Canopy Quality", ylab = "Recruitment",
          frame = FALSE)
  lines(seq(min(can), max(can), length.out=npred), apply(pred.r.can,2,function(x) quantile(x, 0.975)), col="red", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(can), max(can), length.out=npred), apply(pred.r.can,2,function(x) quantile(x, 0.025)), col="red", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(cm0, apply(pred.r.cm,2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(-15, 4.25), col="red",
          xlab = "Channel Modification (scaled)", ylab = "",
          frame = FALSE)
  lines(seq(min(cm), max(cm), length.out=npred), apply(pred.r.cm,2,function(x) quantile(x, 0.975)), col="red", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(cm), max(cm), length.out=npred), apply(pred.r.cm,2,function(x) quantile(x, 0.025)), col="red", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(tq0, apply(pred.r.tq,2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0.15, 0.4), col="black",
          xlab = "TQavg (% Year Above Mean)", ylab = "",
          frame = FALSE)
  lines(seq(min(tq), max(tq), length.out=npred), apply(pred.r.tq,2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(tq), max(tq), length.out=npred), apply(pred.r.tq,2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(rz0, apply(pred.r.rz,2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(5, 10), col="black",
          xlab = "Riparrian zone", ylab = "",
          frame = FALSE)
  lines(seq(min(rz), max(rz), length.out=npred), apply(pred.r.rz,2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(rz), max(rz), length.out=npred), apply(pred.r.rz,2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  
  ## --------------------------------------------------
  # Plot predicted probability a site has adults and NO juveniles (i.e., RELICT; psi*(1-r)) 
  # relationships can be non-linear. 
  
  par(mfrow = c(2,4), mar=c(5,5,5,2), cex.lab = 1.75, cex.axis = 1.5)
  matplot(sub0, apply(pred.psi.sub*(1-pred.r.sub),2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0, 15), col="black",
          xlab = "Substrate Quality", ylab = "Pr(Adults only)", 
          frame = FALSE)
  lines(seq(min(sub), max(sub), length.out=npred), apply(pred.psi.sub*(1-pred.r.sub),2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(sub), max(sub), length.out=npred), apply(pred.psi.sub*(1-pred.r.sub),2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(cover0, apply(pred.psi.cover*(1-pred.r.cover),2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(5, 20), col="black",
          xlab = "Cover Quality", ylab = "",
          frame = FALSE)
  lines(seq(min(cover), max(cover), length.out=npred), apply(pred.psi.cover*(1-pred.r.cover),2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(cover), max(cover), length.out=npred), apply(pred.psi.cover*(1-pred.r.cover),2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(pv0, apply(pred.psi.pv*(1-pred.r.pv),2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0, 10), col="black",
          xlab = "Pool Variety", ylab = "",
          frame = FALSE)
  lines(seq(min(pv), max(pv), length.out=npred), apply(pred.psi.pv*(1-pred.r.pv),2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(pv), max(pv), length.out=npred), apply(pred.psi.pv*(1-pred.r.pv),2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(bs0, apply(pred.psi.bs*(1-pred.r.bs),2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(-8, 6), col="black",
          xlab = "Bank Stability (scaled)", ylab = "",
          frame = FALSE)
  lines(seq(min(bs), max(bs), length.out=npred), apply(pred.psi.bs*(1-pred.r.bs),2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(bs), max(bs), length.out=npred), apply(pred.psi.bs*(1-pred.r.bs),2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(can0,apply(pred.psi.can*(1-pred.r.can),2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0, 10), col="black",
          xlab = "Canopy Quality", ylab = "Pr(Adults only)",
          frame = FALSE)
  lines(seq(min(can), max(can), length.out=npred), apply(pred.psi.can*(1-pred.r.can),2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(can), max(can), length.out=npred), apply(pred.psi.can*(1-pred.r.can),2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(cm0, apply(pred.psi.cm*(1-pred.r.cm),2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(-15, 4.25), col="black",
          xlab = "Channel Modification (scaled)", ylab = "",
          frame = FALSE)
  lines(seq(min(cm), max(cm), length.out=npred), apply(pred.psi.cm*(1-pred.r.cm),2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(cm), max(cm), length.out=npred), apply(pred.psi.cm*(1-pred.r.cm),2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(tq0, apply(pred.psi.tq*(1-pred.r.tq),2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0.15, 0.4), col="black",
          xlab = "TQavg (% Year Above Mean)", ylab = "",
          frame = FALSE)
  lines(seq(min(tq), max(tq), length.out=npred), apply(pred.psi.tq*(1-pred.r.tq),2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(tq), max(tq), length.out=npred), apply(pred.psi.tq*(1-pred.r.tq),2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(rz0, apply(pred.psi.rz*(1-pred.r.rz),2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(5, 10), col="black",
          xlab = "Riparrian zone", ylab = "",
          frame = FALSE)
  lines(seq(min(rz), max(rz), length.out=npred), apply(pred.psi.rz*(1-pred.r.rz),2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(rz), max(rz), length.out=npred), apply(pred.psi.rz*(1-pred.r.rz),2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI

  ## --------------------------------------------------
  # Plot predicted probability a site has adults and juveniles (psi*r) 
  # relationships can be non-linear. 
  
  par(mfrow = c(2,4), mar=c(5,5,5,2), cex.lab = 1.75, cex.axis = 1.5)
  matplot(sub0, apply(pred.psi.sub*pred.r.sub,2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0, 15), col="black",
          xlab = "Substrate Quality", ylab = "Pr(Adults&Juveniles)", 
          frame = FALSE)
  lines(seq(min(sub), max(sub), length.out=npred), apply(pred.psi.sub*pred.r.sub,2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(sub), max(sub), length.out=npred), apply(pred.psi.sub*pred.r.sub,2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(cover0, apply(pred.psi.cover*pred.r.cover,2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(5, 20), col="black",
          xlab = "Cover Quality", ylab = "",
          frame = FALSE)
  lines(seq(min(cover), max(cover), length.out=npred), apply(pred.psi.cover*pred.r.cover,2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(cover), max(cover), length.out=npred), apply(pred.psi.cover*pred.r.cover,2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(pv0, apply(pred.psi.pv*pred.r.pv,2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0, 10), col="black",
          xlab = "Pool Variety", ylab = "",
          frame = FALSE)
  lines(seq(min(pv), max(pv), length.out=npred), apply(pred.psi.pv*pred.r.pv,2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(pv), max(pv), length.out=npred), apply(pred.psi.pv*pred.r.pv,2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(bs0, apply(pred.psi.bs*pred.r.bs,2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(-8, 6), col="black",
          xlab = "Bank Stability (scaled)", ylab = "",
          frame = FALSE)
  lines(seq(min(bs), max(bs), length.out=npred), apply(pred.psi.bs*pred.r.bs,2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(bs), max(bs), length.out=npred), apply(pred.psi.bs*pred.r.bs,2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(can0,apply(pred.psi.can*pred.r.can,2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0, 10), col="black",
          xlab = "Canopy Quality", ylab = "Pr(Adults&Juveniles)",
          frame = FALSE)
  lines(seq(min(can), max(can), length.out=npred), apply(pred.psi.can*pred.r.can,2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(can), max(can), length.out=npred), apply(pred.psi.can*pred.r.can,2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(cm0, apply(pred.psi.cm*pred.r.cm,2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(-15, 4.25), col="black",
          xlab = "Channel Modification (scaled)", ylab = "",
          frame = FALSE)
  lines(seq(min(cm), max(cm), length.out=npred), apply(pred.psi.cm*pred.r.cm,2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(cm), max(cm), length.out=npred), apply(pred.psi.cm*pred.r.cm,2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(tq0, apply(pred.psi.tq*pred.r.tq,2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0.15, 0.4), col="black",
          xlab = "TQavg (% Year Above Mean)", ylab = "",
          frame = FALSE)
  lines(seq(min(tq), max(tq), length.out=npred), apply(pred.psi.tq*pred.r.tq,2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(tq), max(tq), length.out=npred), apply(pred.psi.tq*pred.r.tq,2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
  
  matplot(rz0, apply(pred.psi.rz*pred.r.rz,2,mean), type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(5, 10), col="black",
          xlab = "Riparian zone", ylab = "",
          frame = FALSE)
  lines(seq(min(rz), max(rz), length.out=npred), apply(pred.psi.rz*pred.r.rz,2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(seq(min(rz), max(rz), length.out=npred), apply(pred.psi.rz*pred.r.rz,2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI

## --------------------------------------------------
## Prediction plots - Quadratic of discharge on detection
## --------------------------------------------------  

## !!!IMPORTANT!!!
## Discharge is extremely skewed
  range(discharge); mean(discharge); median(discharge)
  hist(discharge, breaks=seq(0,750,25))
  
# Make predictions across range of observed values
  npred <- 500
  pred.p2.d<-pred.p32.d<-pred.p33.d<-matrix(NA, nrow(outmat), npred)
  
  #Discharge
  d0 <- seq(min(discharge), max(discharge), length.out=npred)
  dP <- standardize2match(d0, discharge)
  
  #Discharge^2
  d02<-d0^2
  dP2<-standardize2match(d02, discharge2)  

## As we did above, write out the JAGS model to make predictions during each iteration
  for(i in 1:nrow(outmat)){
    ## predict probabilities for state 2 (adults only)
    pred.p2.d[i,] <- plogis(qlogis(outmat[i,"mean.p2"]) + outmat[i,"beta.lp2[1]"] * dP + outmat[i,"beta.lp2[2]"] * dP2 )
    
    ## predict probabilities for 
    ##    state 3 (adults and juvs) - observe only adults
    ##    state 3 (adults and juvs) - observe adults and juvs
    
    ## multinomial logit
    m.pred.p32.d <- outmat[i,"alpha.lp32"] + outmat[i,"beta.lp32[1]"] * dP + outmat[i,"beta.lp32[2]"] * dP2
    m.pred.p33.d <- outmat[i,"alpha.lp33"] + outmat[i,"beta.lp33[1]"] * dP + outmat[i,"beta.lp33[2]"] * dP2
    
    pred.p32.d[i,] <- exp(m.pred.p32.d) / (1 + exp(m.pred.p32.d) + exp(m.pred.p33.d))
    pred.p33.d[i,] <- exp(m.pred.p33.d) / (1 + exp(m.pred.p32.d) + exp(m.pred.p33.d))
    
  }    
  
  
## Plots
par(mfrow = c(2,2), mar=c(5,5,5,2), cex.lab = 1.75, cex.axis = 1.5)

## Detection of adults if site is only occupied by adults
  matplot(d0,apply(pred.p2.d, 2, mean), type = "l", lwd = 2, ylim = c(0,1), xlim = c(0, 300),
        xlab = "Discharge m3/s", ylab = "Adults at Relict",
        frame = FALSE)
  lines(d0, apply(pred.p2.d,2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(d0, apply(pred.p2.d,2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
  
## Detection of adults if site is occupied by adults and juvs
  matplot(d0, apply(pred.p32.d, 2, mean), typ="l", lwd = 2, ylim = c(0,1), xlim = c(0,300), 
        xlab="Discharge m3/s", ylab="Adults at Recruit.",
        frame = FALSE)
  lines(d0, apply(pred.p32.d,2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(d0, apply(pred.p32.d,2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
  
## Detection of adults and juvs if site is occupied by adults and juvs
  matplot(d0, apply(pred.p33.d, 2, mean), typ="l", lwd = 2, ylim = c(0,1), xlim = c(0, 300),
        xlab="Discharge m3/s", ylab="Juv. at Recruit",
        frame = FALSE)
  lines(d0, apply(pred.p33.d,2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(d0, apply(pred.p33.d,2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
 
## Prob of not detecting adults or and adults and juvs if site is occupied by adults and juvs
  matplot(d0, apply(1-(pred.p32.d+pred.p33.d), 2, mean), typ="l", lwd = 2, ylim = c(0,1), xlim = c(0, 300),
        xlab="Discharge m3/s", ylab="Failure at Recruit.",
        frame = FALSE)
  lines(d0, apply(1-(pred.p32.d+pred.p33.d),2,function(x) quantile(x, 0.975)), col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
  lines(d0, apply(1-(pred.p32.d+pred.p33.d),2,function(x) quantile(x, 0.025)), col="black", pch=22, lty=2, lwd=2) #lower 95% CRI


