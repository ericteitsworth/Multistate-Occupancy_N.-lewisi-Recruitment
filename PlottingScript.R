#Teitsworth et al. 

## Plotting Multistate Model Results ##


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


### Plot predictions for supported and partially supported covariate effects in psi and r ###
# psi -> substrate, cover, pool variety, bank stability, canopy, channel mod, TQavg
# r -> cover, pool variety, riparian, TQavg

#grab posterior means
tmp<- model$mean

npred<-500

pred.psi.sub<-pred.psi.cover<-pred.psi.pv<-pred.psi.bs<-pred.psi.can<-pred.psi.cm<-pred.psi.tq<-pred.r.cover<-pred.r.tq<-pred.r.pv<-pred.r.rz<-array(NA, dim=c(npred,1))

sub0 <- seq(min(sub), max(sub), length.out=npred)
cover0 <- seq(min(cover), max(cover), length.out=npred)
cm0 <- seq(min(cm), max(cm), length.out=npred)
pv0 <- seq(min(pv), max(pv), length.out=npred)
rz0 <- seq(min(rz), max(rz), length.out=npred)
tq0 <- seq(min(tq), max(tq), length.out=npred)
bs0 <- seq(min(bs), max(bs), length.out=npred)
can0 <- seq(min(can), max(can), length.out=npred)
require(AHMbook)
subP <- standardize2match(sub0, sub)
coverP <- standardize2match(cover0, cover)
cmP <- standardize2match(cm0, cm)
pvP <- standardize2match(pv0, pv)
rzP <- standardize2match(rz0, rz)
tqP <- standardize2match(tq0, tq)
bsP <- standardize2match(bs0, bs)
canP <- standardize2match(can0, can)

# Predict at the means of the other covariate (i.e., 0)
pred.psi.sub <- plogis(tmp$mean.psi + tmp$beta.lpsi[1] * subP)
pred.psi.cover <- plogis(tmp$mean.psi + tmp$beta.lpsi[2] * coverP)
pred.psi.cm <- plogis(tmp$mean.psi + tmp$beta.lpsi[4] * cmP)
pred.psi.pv <- plogis(tmp$mean.psi + tmp$beta.lpsi[7] * pvP)
pred.psi.rz <- plogis(tmp$mean.psi + tmp$beta.lpsi[6] * rzP)
pred.psi.tq <- plogis(tmp$mean.psi + tmp$beta.lpsi[3] * tqP)
pred.psi.bs <- plogis(tmp$mean.psi + tmp$beta.lpsi[8] * bsP)
pred.psi.can <- plogis(tmp$mean.psi + tmp$beta.lpsi[5] * canP)
pred.r.sub <- plogis(tmp$mean.r + tmp$beta.lr[1] * subP)
pred.r.cover <- plogis(tmp$mean.r + tmp$beta.lr[2] * coverP)
pred.r.cm <- plogis(tmp$mean.r + tmp$beta.lr[4] * cmP)
pred.r.pv <- plogis(tmp$mean.r + tmp$beta.lr[7] * pvP)
pred.r.rz <- plogis(tmp$mean.r + tmp$beta.lr[6] * rzP)
pred.r.tq <- plogis(tmp$mean.r + tmp$beta.lr[3] * tqP)
pred.r.bs <- plogis(tmp$mean.r + tmp$beta.lr[8] * bsP)
pred.r.can <- plogis(tmp$mean.r + tmp$beta.lr[5] * canP)

######## Upper 95% CRI
# Substrate
orig.pred.sub.up <- seq(min(sub), max(sub), length.out=npred)
pred.sub.up <- (orig.pred.sub.up - mean(sub)) / sd(sub)
pred.psi.sub.up <- array(NA, dim = c(length(pred.sub.up), 1))
for (i in 1:length(pred.sub.up)){
  pred.psi.sub.up[i,] <- plogis((model$summary[1,1] + (model$summary[2,2]*1.96)) + model$summary[2,1] * pred.sub.up[i])
}

# Cover
orig.pred.cover.up <- seq(min(cover), max(cover), length.out=npred)
pred.cover.up <- (orig.pred.cover.up - mean(cover)) / sd(cover)
pred.psi.cover.up <- array(NA, dim = c(length(pred.cover.up), 1))
for (i in 1:length(pred.cover.up)){
  pred.psi.cover.up[i,] <- plogis((model$summary[1,1] + (model$summary[3,2]*1.96)) + model$summary[3,1] * pred.cover.up[i])
}
pred.r.cover.up <- array(NA, dim = c(length(pred.cover.up), 1))
for (i in 1:length(pred.cover.up)){
  pred.r.cover.up[i,] <- plogis((model$summary[10,1] + (model$summary[12,2]*1.96)) + model$summary[12,1] * pred.cover.up[i])
}

# Channel
orig.pred.cm.up <- seq(min(cm), max(cm), length.out=npred)
pred.cm.up <- (orig.pred.cm.up - mean(cm)) / sd(cm)
pred.psi.cm.up <- array(NA, dim = c(length(pred.cm.up), 1))
for (i in 1:length(pred.cm.up)){
  pred.psi.cm.up[i,] <- plogis((model$summary[1,1] + (model$summary[5,2]*1.96)) + model$summary[5,1] * pred.cm.up[i])
}

# Pool Variety
orig.pred.pv.up <- seq(min(pv), max(pv), length.out=npred)
pred.pv.up <- (orig.pred.pv.up - mean(pv)) / sd(pv)
pred.psi.pv.up <- array(NA, dim = c(length(pred.pv.up), 1))
for (i in 1:length(pred.pv.up)){
  pred.psi.pv.up[i,] <- plogis((model$summary[1,1] + (model$summary[8,2]*1.96)) + model$summary[8,1] * pred.pv.up[i])
}
pred.r.pv.up <- array(NA, dim = c(length(pred.pv.up), 1))
for (i in 1:length(pred.pv.up)){
  pred.r.pv.up[i,] <- plogis((model$summary[10,1] + (model$summary[17,2]*1.96)) + model$summary[17,1] * pred.pv.up[i])
}

# Riparian Zone
orig.pred.rz.up <- seq(min(rz), max(rz), length.out=npred)
pred.rz.up <- (orig.pred.rz.up - mean(rz)) / sd(rz)
pred.r.rz.up <- array(NA, dim = c(length(pred.rz.up), 1))
for (i in 1:length(pred.rz.up)){
  pred.r.rz.up[i,] <- plogis((model$summary[10,1] + (model$summary[16,2]*1.96)) + model$summary[16,1] * pred.rz.up[i])
}

# TQavg
orig.pred.tq.up <- seq(min(tq), max(tq), length.out=npred)
pred.tq.up <- (orig.pred.tq.up - mean(tq)) / sd(tq)
pred.psi.tq.up <- array(NA, dim = c(length(pred.tq.up), 1))
for (i in 1:length(pred.tq.up)){
  pred.psi.tq.up[i,] <- plogis((model$summary[1,1] + (model$summary[4,2]*1.96)) + model$summary[4,1] * pred.tq.up[i])
}
pred.r.tq.up <- array(NA, dim = c(length(pred.tq.up), 1))
for (i in 1:length(pred.tq.up)){
  pred.r.tq.up[i,] <- plogis((model$summary[10,1] + (model$summary[13,2]*1.96)) + model$summary[13,1] * pred.tq.up[i])
}

# Bank Stability
orig.pred.bs.up <- seq(min(bs), max(bs), length.out=npred)
pred.bs.up <- (orig.pred.bs.up - mean(bs)) / sd(bs)
pred.psi.bs.up <- array(NA, dim = c(length(pred.bs.up), 1))
for (i in 1:length(pred.bs.up)){
  pred.psi.bs.up[i,] <- plogis((model$summary[1,1] + (model$summary[9,2]*1.96)) + model$summary[9,1] * pred.bs.up[i])
}

# Canopy Cover
orig.pred.can.up <- seq(min(can), max(can), length.out=npred)
pred.can.up <- (orig.pred.can.up - mean(can)) / sd(can)
pred.psi.can.up <- array(NA, dim = c(length(pred.can.up), 1))
for (i in 1:length(pred.can.up)){
  pred.psi.can.up[i,] <- plogis((model$summary[1,1] + (model$summary[6,2]*1.96)) + model$summary[6,1] * pred.can.up[i])
}


########## Lower 95% CRI
# Substrate
orig.pred.sub.down <- seq(min(sub), max(sub), length.out=npred)
pred.sub.down <- (orig.pred.sub.down - mean(sub)) / sd(sub)
pred.psi.sub.down <- array(NA, dim = c(length(pred.sub.down), 1))
for (i in 1:length(pred.sub.down)){
  pred.psi.sub.down[i,] <- plogis((model$summary[1,1] - (model$summary[2,2]*1.96)) + model$summary[2,1] * pred.sub.down[i])
}

# Cover
orig.pred.cover.down <- seq(min(cover), max(cover), length.out=npred)
pred.cover.down <- (orig.pred.cover.down - mean(cover)) / sd(cover)
pred.psi.cover.down <- array(NA, dim = c(length(pred.cover.down), 1))
for (i in 1:length(pred.cover.down)){
  pred.psi.cover.down[i,] <- plogis((model$summary[1,1] - (model$summary[3,2]*1.96)) + model$summary[3,1] * pred.cover.down[i])
}
pred.r.cover.down <- array(NA, dim = c(length(pred.cover.down), 1))
for (i in 1:length(pred.cover.down)){
  pred.r.cover.down[i,] <- plogis((model$summary[10,1] - (model$summary[12,2]*1.96)) + model$summary[12,1] * pred.cover.down[i])
}

# Channel
orig.pred.cm.down <- seq(min(cm), max(cm), length.out=npred)
pred.cm.down <- (orig.pred.cm.down - mean(cm)) / sd(cm)
pred.psi.cm.down <- array(NA, dim = c(length(pred.cm.down), 1))
for (i in 1:length(pred.cm.down)){
  pred.psi.cm.down[i,] <- plogis((model$summary[1,1] - (model$summary[5,2]*1.96)) + model$summary[5,1] * pred.cm.down[i])
}

# Pool Variety
orig.pred.pv.down <- seq(min(pv), max(pv), length.out=npred)
pred.pv.down <- (orig.pred.pv.down - mean(pv)) / sd(pv)
pred.psi.pv.down <- array(NA, dim = c(length(pred.pv.down), 1))
for (i in 1:length(pred.pv.down)){
  pred.psi.pv.down[i,] <- plogis((model$summary[1,1] - (model$summary[8,2]*1.96)) + model$summary[8,1] * pred.pv.down[i])
}
pred.r.pv.down <- array(NA, dim = c(length(pred.pv.down), 1))
for (i in 1:length(pred.pv.down)){
  pred.r.pv.down[i,] <- plogis((model$summary[10,1] - (model$summary[17,2]*1.96)) + model$summary[17,1] * pred.pv.down[i])
}

# Riparian Zone
orig.pred.rz.down <- seq(min(rz), max(rz), length.out=npred)
pred.rz.down <- (orig.pred.rz.down - mean(rz)) / sd(rz)
pred.r.rz.down <- array(NA, dim = c(length(pred.rz.down), 1))
for (i in 1:length(pred.rz.down)){
  pred.r.rz.down[i,] <- plogis((model$summary[10,1] - (model$summary[16,2]*1.96)) + model$summary[16,1] * pred.rz.down[i])
}

# TQavg
orig.pred.tq.down <- seq(min(tq), max(tq), length.out=npred)
pred.tq.down <- (orig.pred.tq.down - mean(tq)) / sd(tq)
pred.psi.tq.down <- array(NA, dim = c(length(pred.tq.down), 1))
for (i in 1:length(pred.tq.down)){
  pred.psi.tq.down[i,] <- plogis((model$summary[1,1] - (model$summary[4,2]*1.96)) + model$summary[4,1] * pred.tq.down[i])
}
pred.r.tq.down <- array(NA, dim = c(length(pred.tq.down), 1))
for (i in 1:length(pred.tq.down)){
  pred.r.tq.down[i,] <- plogis((model$summary[10,1] - (model$summary[13,2]*1.96)) + model$summary[13,1] * pred.tq.down[i])
}

# Bank Stability
orig.pred.bs.down <- seq(min(bs), max(bs), length.out=npred)
pred.bs.down <- (orig.pred.bs.down - mean(bs)) / sd(bs)
pred.psi.bs.down <- array(NA, dim = c(length(pred.bs.down), 1))
for (i in 1:length(pred.bs.down)){
  pred.psi.bs.down[i,] <- plogis((model$summary[1,1] - (model$summary[9,2]*1.96)) + model$summary[9,1] * pred.bs.down[i])
}

# Canopy Cover
orig.pred.can.down <- seq(min(can), max(can), length.out=npred)
pred.can.down <- (orig.pred.can.down - mean(can)) / sd(can)
pred.psi.can.down <- array(NA, dim = c(length(pred.can.down), 1))
for (i in 1:length(pred.can.down)){
  pred.psi.can.down[i,] <- plogis((model$summary[1,1] - (model$summary[6,2]*1.96)) + model$summary[6,1] * pred.can.down[i])
}


#### Plot supported and weakly supported covariate effects of psi and r ####
#Plot psi effects
par(mfrow = c(2,4), mar=c(5,5,5,2), cex.lab = 1.75, cex.axis = 1.5)
matplot(sub0, pred.psi.sub, type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0, 15), col="red",
        xlab = "Substrate Quality", ylab = "Occupancy", 
        #main = "Substrate Quality effect on Occupancy", 
        frame = FALSE)
lines(seq(min(sub), max(sub), length.out=npred), pred.psi.sub.up, col="red", pch=22, lty=2, lwd=2) #upper 95% CRI
lines(seq(min(sub), max(sub), length.out=npred), pred.psi.sub.down, col="red", pch=22, lty=2, lwd=2) #lower 95% CRI
matplot(cover0, pred.psi.cover, type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(5, 20), col="red",
        xlab = "Cover Quality", ylab = "",
        #main = "Cover Quality effect on Occupancy", 
        frame = FALSE)
lines(seq(min(cover), max(cover), length.out=npred), pred.psi.cover.up, col="red", pch=22, lty=2, lwd=2) #upper 95% CRI
lines(seq(min(cover), max(cover), length.out=npred), pred.psi.cover.down, col="red", pch=22, lty=2, lwd=2) #lower 95% CRI
matplot(pv0, pred.psi.pv, type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0, 10), col="red",
        xlab = "Pool Variety", ylab = "",
        # main = "Pool Variety effect on Occupancy", 
        frame = FALSE)
lines(seq(min(pv), max(pv), length.out=npred), pred.psi.pv.up, col="red", pch=22, lty=2, lwd=2) #upper 95% CRI
lines(seq(min(pv), max(pv), length.out=npred), pred.psi.pv.down, col="red", pch=22, lty=2, lwd=2) #lower 95% CRI
matplot(bs0, pred.psi.bs, type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0, 6),
        xlab = "Bank Stability", ylab = "",
        # main = "Bank Stability effect on Occupancy", 
        frame = FALSE)
lines(seq(min(bs), max(bs), length.out=npred), pred.psi.bs.up, col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
lines(seq(min(bs), max(bs), length.out=npred), pred.psi.bs.down, col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
matplot(can0, pred.psi.can, type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0, 10),
        xlab = "Canopy Quality", ylab = "Occupancy",
        # main = "Canopy Quality effect on Occupancy", 
        frame = FALSE)
lines(seq(min(can), max(can), length.out=npred), pred.psi.can.up, col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
lines(seq(min(can), max(can), length.out=npred), pred.psi.can.down, col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
matplot(cm0, pred.psi.cm, type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(-15, 4.25), col="red",
        xlab = "Channel Modification (scaled)", ylab = "",
        # main = "Channel Modification effect on Occupancy", 
        frame = FALSE)
lines(seq(min(cm), max(cm), length.out=npred), pred.psi.cm.up, col="red", pch=22, lty=2, lwd=2) #upper 95% CRI
lines(seq(min(cm), max(cm), length.out=npred), pred.psi.cm.down, col="red", pch=22, lty=2, lwd=2) #lower 95% CRI
matplot(tq0, pred.psi.tq, type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0.15, 0.4),
        xlab = "TQavg (% Year Above Mean)", ylab = "",
        #main = "Stream Flashiness effect on Occupancy", 
        frame = FALSE)
lines(seq(min(tq), max(tq), length.out=npred), pred.psi.tq.up, col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
lines(seq(min(tq), max(tq), length.out=npred), pred.psi.tq.down, col="black", pch=22, lty=2, lwd=2) #lower 95% CRI

#Plot r effects
par(mfrow = c(1,4), mar=c(5,5,5,4), cex.lab = 1.75, cex.axis = 1.5)  
matplot(cover0, pred.r.cover, type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(5, 20),
        xlab = "Cover Quality", ylab = "Recruitment",
        #main = "Cover Availability effect on Recruitment", 
        frame = FALSE)
lines(seq(min(cover), max(cover), length.out=npred), pred.r.cover.up, col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
lines(seq(min(cover), max(cover), length.out=npred), pred.r.cover.down, col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
matplot(pv0, pred.r.pv, type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0,10),
        xlab = "Pool Variety", ylab = "",
        #main = "Pool Variety effect on Recruitment", 
        frame = FALSE)
lines(seq(min(pv), max(pv), length.out=npred), pred.r.pv.up, col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
lines(seq(min(pv), max(pv), length.out=npred), pred.r.pv.down, col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
matplot(rz0, pred.r.rz, type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(5,10),
        xlab = "Riparian Zone Quality", ylab = "",
        # main = "Riparian Quality effect on Recruitment", 
        frame = FALSE)
lines(seq(min(rz), max(rz), length.out=npred), pred.r.rz.up, col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
lines(seq(min(rz), max(rz), length.out=npred), pred.r.rz.down, col="black", pch=22, lty=2, lwd=2) #lower 95% CRI  
matplot(tq0, pred.r.tq, type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0.15, 0.4),
        xlab = "TQavg (% Year Above Mean)", ylab = "",
        #main = "Stream Flashiness effect on Recruitment",
        frame = FALSE)
lines(seq(min(tq), max(tq), length.out=npred), pred.r.tq.up, col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
lines(seq(min(tq), max(tq), length.out=npred), pred.r.tq.down, col="black", pch=22, lty=2, lwd=2) #lower 95% CRI



###############################################################################################

#Plot discharge effects on detection
npred <- 500
pred.p2.d<-pred.p32.d<-pred.p33.d<-array(NA, dim=c(npred,1))

#Discharge
d0 <- seq(min(discharge), max(discharge), length.out=npred)
dP <- standardize2match(d0, discharge)
pred.p2.d <- plogis(tmp$mean.p2 + tmp$beta.lp2[1] * dP)
pred.p32.d <- plogis(tmp$alpha.lp32 + tmp$beta.lp32[1] * dP)
pred.p33.d <- plogis(tmp$alpha.lp33 + tmp$beta.lp33[1] * dP)

#Plot quadratic effects of discharge on detection p33 and p32
pred.p33.d2 <- pred.p32.d2<- array(NA, dim=c(npred,1))
d20 <- seq(min(d2.scaled), max(d2.scaled), length.out=npred)
d2P <- standardize2match(d20, d2.scaled)
pred.p33.d2 <- plogis(tmp$alpha.lp33 + tmp$beta.lp33[2] * d2P)
pred.p32.d2 <- plogis(tmp$alpha.lp32 + tmp$beta.lp32[2] * d2P)

#Bait Age
pred.p32.b<-pred.p33.b<- array(NA, dim=c(npred,1))
b0 <- seq(min(bait), max(bait), length.out=npred)
bP <- standardize2match(b0, bait)
pred.p32.b <- plogis(tmp$alpha.lp32 + tmp$beta.lp32[3] * bP)
pred.p33.b <- plogis(tmp$alpha.lp33 + tmp$beta.lp33[3] * bP)

######### Upper 95% CRI
# Discharge on p2
orig.pred.dp2.up <- seq(min(discharge), max(discharge), length.out=npred)
pred.dp2.up <- (orig.pred.dp2.up - mean(discharge)) / sd(discharge)
pred.p2.up <- array(NA, dim = c(length(pred.dp2.up), 1))
for (i in 1:length(pred.dp2.up)){
  pred.p2.up[i,] <- plogis((model$summary[19,1] + (model$summary[22,2]*1.96)) + model$summary[22,1] * pred.dp2.up[i])
}
# Discharge on p32
orig.pred.dp32.up <- seq(min(discharge), max(discharge), length.out=npred)
pred.dp32.up <- (orig.pred.dp32.up - mean(discharge)) / sd(discharge)
pred.p32.up <- array(NA, dim = c(length(pred.dp32.up), 1))
for (i in 1:length(pred.dp32.up)){
  pred.p32.up[i,] <- plogis((model$summary[20,1] + (model$summary[25,2]*1.96)) + model$summary[25,1] * pred.dp32.up[i])
}
# Discharge on p33
orig.pred.dp33.up <- seq(min(discharge), max(discharge), length.out=npred)
pred.dp33.up <- (orig.pred.dp33.up - mean(discharge)) / sd(discharge)
pred.p33.up <- array(NA, dim = c(length(pred.dp33.up), 1))
for (i in 1:length(pred.dp33.up)){
  pred.p33.up[i,] <- plogis((model$summary[21,1] + (model$summary[28,2]*1.96)) + model$summary[28,1] * pred.dp33.up[i])
}
# Discharge^2 on p33
orig.pred.d2p33.up <- seq(min(d2.scaled), max(d2.scaled), length.out=npred)
pred.d2p33.up <- (orig.pred.d2p33.up - mean(d2.scaled)) / sd(d2.scaled)
pred.2p33.up <- array(NA, dim = c(length(pred.d2p33.up), 1))
for (i in 1:length(pred.d2p33.up)){
  pred.2p33.up[i,] <- plogis((model$summary[21,1] + (model$summary[29,2]*1.96)) + model$summary[29,1] * pred.d2p33.up[i])
}
# Discharge^2 on p32
orig.pred.d2p32.up <- seq(min(d2.scaled), max(d2.scaled), length.out=npred)
pred.d2p32.up <- (orig.pred.d2p32.up - mean(d2.scaled)) / sd(d2.scaled)
pred.2p32.up <- array(NA, dim = c(length(pred.d2p32.up), 1))
for (i in 1:length(pred.d2p32.up)){
  pred.2p32.up[i,] <- plogis((model$summary[20,1] + (model$summary[26,2]*1.96)) + model$summary[26,1] * pred.d2p32.up[i])
}
# Bait Age on p33
orig.pred.bp33.up <- seq(min(bait), max(bait), length.out=npred)
pred.bp33.up <- (orig.pred.bp33.up - mean(bait)) / sd(bait)
pred.bait.p33.up <- array(NA, dim = c(length(pred.bp33.up), 1))
for (i in 1:length(pred.bp33.up)){
  pred.bait.p33.up[i,] <- plogis((model$summary[21,1] + (model$summary[30,2]*1.96)) + model$summary[30,1] * pred.bp33.up[i])
}
# Bait Age on p32
orig.pred.bp32.up <- seq(min(bait), max(bait), length.out=npred)
pred.bp32.up <- (orig.pred.bp32.up - mean(bait)) / sd(bait)
pred.bait.p32.up <- array(NA, dim = c(length(pred.bp32.up), 1))
for (i in 1:length(pred.bp32.up)){
  pred.bait.p32.up[i,] <- plogis((model$summary[20,1] + (model$summary[27,2]*1.96)) + model$summary[27,1] * pred.bp32.up[i])
}


### Lower 95% CRI
# Discharge on p2
orig.pred.dp2.down <- seq(min(discharge), max(discharge), length.out=npred)
pred.dp2.down <- (orig.pred.dp2.down - mean(discharge)) / sd(discharge)
pred.p2.down <- array(NA, dim = c(length(pred.dp2.down), 1))
for (i in 1:length(pred.dp2.down)){
  pred.p2.down[i,] <- plogis((model$summary[19,1] - (model$summary[22,2]*1.96)) + model$summary[22,1] * pred.dp2.down[i])
}
# Discharge on p32
orig.pred.dp32.down <- seq(min(discharge), max(discharge), length.out=npred)
pred.dp32.down <- (orig.pred.dp32.down - mean(discharge)) / sd(discharge)
pred.p32.down <- array(NA, dim = c(length(pred.dp32.down), 1))
for (i in 1:length(pred.dp32.down)){
  pred.p32.down[i,] <- plogis((model$summary[20,1] - (model$summary[25,2]*1.96)) + model$summary[25,1] * pred.dp32.down[i])
}
# Discharge on p33
orig.pred.dp33.down <- seq(min(discharge), max(discharge), length.out=npred)
pred.dp33.down <- (orig.pred.dp33.down - mean(discharge)) / sd(discharge)
pred.p33.down <- array(NA, dim = c(length(pred.dp33.down), 1))
for (i in 1:length(pred.dp33.down)){
  pred.p33.down[i,] <- plogis((model$summary[21,1] - (model$summary[28,2]*1.96)) + model$summary[28,1] * pred.dp33.down[i])
}
# Discharge^2 on p33
orig.pred.d2p33.down <- seq(min(d2.scaled), max(d2.scaled), length.out=npred)
pred.d2p33.down <- (orig.pred.d2p33.down - mean(d2.scaled)) / sd(d2.scaled)
pred.2p33.down <- array(NA, dim = c(length(pred.d2p33.down), 1))
for (i in 1:length(pred.d2p33.down)){
  pred.2p33.down[i,] <- plogis((model$summary[21,1] - (model$summary[29,2]*1.96)) + model$summary[29,1] * pred.d2p33.down[i])
}
# Discharge^2 on p32
orig.pred.d2p32.down <- seq(min(d2.scaled), max(d2.scaled), length.out=npred)
pred.d2p32.down <- (orig.pred.d2p32.down - mean(d2.scaled)) / sd(d2.scaled)
pred.2p32.down <- array(NA, dim = c(length(pred.d2p32.down), 1))
for (i in 1:length(pred.d2p32.down)){
  pred.2p32.down[i,] <- plogis((model$summary[20,1] - (model$summary[26,2]*1.96)) + model$summary[26,1] * pred.d2p32.down[i])
}
# Bait Age on p33
orig.pred.bp33.down <- seq(min(bait), max(bait), length.out=npred)
pred.bp33.down <- (orig.pred.bp33.down - mean(bait)) / sd(bait)
pred.bait.p33.down <- array(NA, dim = c(length(pred.bp33.down), 1))
for (i in 1:length(pred.bp33.down)){
  pred.bait.p33.down[i,] <- plogis((model$summary[21,1] - (model$summary[30,2]*1.96)) + model$summary[30,1] * pred.bp33.down[i])
}
# Bait Age on p32
orig.pred.bp32.down <- seq(min(bait), max(bait), length.out=npred)
pred.bp32.down <- (orig.pred.bp32.down - mean(bait)) / sd(bait)
pred.bait.p32.down <- array(NA, dim = c(length(pred.bp32.down), 1))
for (i in 1:length(pred.bp32.down)){
  pred.bait.p32.down[i,] <- plogis((model$summary[20,1] - (model$summary[27,2]*1.96)) + model$summary[27,1] * pred.bp32.down[i])
}


#### Plot supported and weakly supported covariate effects of detection ####
par(mfrow = c(2,2), mar = c(5,5,5,4), cex.lab = 1.75, cex.axis = 1.5)
matplot(d0, pred.p32.d, type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0, 775), col="red",
        xlab = "Discharge (m3/s)", ylab = "Adults in r",
        #main = "Discharge effect on Detecting Adults Only", 
        frame = FALSE)
lines(seq(min(discharge), max(discharge), length.out=npred), pred.p32.up, col="red", pch=22, lty=2, lwd=2) #upper 95% CRI
lines(seq(min(discharge), max(discharge), length.out=npred), pred.p32.down, col="red", pch=22, lty=2, lwd=2) #lower 95% CRI
matplot(d0, pred.p33.d, type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0, 775),
        xlab = "Discharge (m3/s)", ylab = "Juveniles in r",
        #main = "Discharge effect on Detecting Juveniles", 
        frame = FALSE)
lines(seq(min(discharge), max(discharge), length.out=npred), pred.p33.up, col="black", pch=22, lty=2, lwd=2) #upper 95% CRI
lines(seq(min(discharge), max(discharge), length.out=npred), pred.p33.down, col="black", pch=22, lty=2, lwd=2) #lower 95% CRI
matplot(d20, pred.p33.d2, type = "l", lty = 1, lwd = 3, ylim = c(0,1), xlim = c(0, 3), col="red",
        xlab = "Discharge^2", ylab = "",
        #main = "Discharge^2 effect on Detecting Juveniles", 
        frame = FALSE)
lines(seq(min(d2.scaled), max(d2.scaled), length.out=npred), pred.2p33.up, col="red", pch=22, lty=2, lwd=2) #upper 95% CRI
lines(seq(min(d2.scaled), max(d2.scaled), length.out=npred), pred.2p33.down, col="red", pch=22, lty=2, lwd=2) #lower 95% CRI
