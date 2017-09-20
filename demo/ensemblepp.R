### R code from vignette source 'Ch11.Rnw'

###################################################
### code chunk number 1: setup
###################################################
set.seed(100)


###################################################
### code chunk number 2: clogispit
###################################################
clogispit <- function(obs, location, scale, left = 0) {
  pit <- pclogis(obs, location, scale, left = left)
  pit[obs <= 0] <- runif(sum(obs <= 0), 0, pit[obs <= 0])
  return(pit)
}


###################################################
### code chunk number 3: mvrank1
###################################################
mv.rank <- function(x) {
  d <- dim(x)
  x.prerank <- numeric(d[2])
  for(i in 1:d[2]) {
    x.prerank[i] <- sum(apply(x<=x[,i],2,all))
  }
  x.rank <- rank(x.prerank,ties="random")
  return(x.rank)
}


###################################################
### code chunk number 4: install (eval = FALSE)
###################################################
## install.packages(c("ensemblepp", "ensembleBMA", "crch", "gamlss", 
## "ensembleMOS",  "SpecsVerification", "scoringRules", "glmx", "ordinal", 
## "pROC", "mvtnorm"))


###################################################
### code chunk number 5: data-temp
###################################################
data("temp", package = "ensemblepp")


###################################################
### code chunk number 6: data-temp2
###################################################
dim(temp)
names(temp)


###################################################
### code chunk number 7: data-temp3
###################################################
temp$date <- as.Date(rownames(temp))
range(temp$date)


###################################################
### code chunk number 8: data-temp4
###################################################
temp <- temp[format(temp$date, "%m") %in% c("12", "01", "02"),]


###################################################
### code chunk number 9: data-temp5
###################################################
temp$ensmean <- apply(temp[,2:12], 1, mean)
temp$enssd   <- apply(temp[,2:12], 1, sd)


###################################################
### code chunk number 10: data-temp6
###################################################
temptrain <-  temp[temp$date <  "2010-03-01",]
temptest  <-  temp[temp$date >= "2010-03-01",]


###################################################
### code chunk number 11: data-temp7
###################################################
dim(temptrain)
dim(temptest)
names(temp)


###################################################
### code chunk number 12: tempscatter
###################################################
plot(temp ~ ensmean, data = temptrain)
abline(0, 1, lty = 2)


###################################################
### code chunk number 13: regressionline
###################################################
MOS <- lm(temp ~ ensmean, data = temptrain)
abline(MOS)


###################################################
### code chunk number 14: fig_tempscatter
###################################################
plot(temp ~ ensmean, data = temptrain)
abline(0, 1, lty = 2)
MOS <- lm(temp ~ ensmean, data = temptrain)
abline(MOS)


###################################################
### code chunk number 15: Ch11.Rnw:266-267
###################################################
summary(MOS)


###################################################
### code chunk number 16: MOSfc
###################################################
fcMOS <- predict(MOS, newdata = temptest)


###################################################
### code chunk number 17: MOSfc2
###################################################
plot(fcMOS[1:20], type = "l", lty = 2, ylab = "2m temperature", 
  xlab = "date", xaxt = "n", ylim = c(-35, 10))
axis(1, at = seq(1,20,6), temptest$date[seq(1, 20, 6)])

lines(temptest$temp[1:20])
lines(temptest$ensmean[1:20], lty = 3)


###################################################
### code chunk number 18: fig_MOSfc
###################################################
plot(fcMOS[1:20], type = "l", lty = 2, ylab = "2m temperature", 
  xlab = "date", xaxt = "n", ylim = c(-35, 10))
axis(1, at = seq(1,20,6), temptest$date[seq(1, 20, 6)])

lines(temptest$temp[1:20])
lines(temptest$ensmean[1:20], lty = 3)


###################################################
### code chunk number 19: MOS-error
###################################################
rbind(raw = c(BIAS =  mean(temptest$ensmean-temptest$temp),
              MAE = mean(abs(temptest$ensmean - temptest$temp)),
              RMSE = sqrt(mean((temptest$ensmean - temptest$temp)^2))),
      MOS = c(BIAS = mean(fcMOS-temptest$temp),
              MAE = mean(abs(fcMOS - temptest$temp)),
              RMSE = sqrt(mean((fcMOS - temptest$temp)^2))))


###################################################
### code chunk number 20: data-temp8
###################################################
library("ensembleBMA")
temptrain_eD <- ensembleData(forecasts = temptrain[,2:12],
  dates = temptrain$date, observations = temptrain$temp,
  forecastHour = 24, initializationTime = "00", exchangeable = rep(1, 11))
temptest_eD <- ensembleData(forecasts = temptest[,2:12],
  dates = temptest$date, observations = temptest$temp,
  forecastHour = 24, initializationTime = "00", exchangeable = rep(1, 11))


###################################################
### code chunk number 21: rankhist
###################################################
  rank <- apply(temptrain[, 1:12], 1, rank)[1,]
  hist(rank, breaks = 0:12 + 0.5, main = "Verification Rank Histogram")


###################################################
### code chunk number 22: spreadskill
###################################################
sdcat <- cut(temptrain$enssd, 
  breaks = quantile(temptrain$enssd, seq(0, 1, 0.2)))
boxplot(abs(residuals(MOS)) ~ sdcat, ylab = "absolute error", 
  xlab = "ensemble standard deviation", main = "Spread-Skill")


###################################################
### code chunk number 23: fig_ensemble
###################################################
par(mfrow = c(1,2))
  rank <- apply(temptrain[, 1:12], 1, rank)[1,]
  hist(rank, breaks = 0:12 + 0.5, main = "Verification Rank Histogram")
sdcat <- cut(temptrain$enssd, 
  breaks = quantile(temptrain$enssd, seq(0, 1, 0.2)))
boxplot(abs(residuals(MOS)) ~ sdcat, ylab = "absolute error", 
  xlab = "ensemble standard deviation", main = "Spread-Skill")


###################################################
### code chunk number 24: NGR-fit
###################################################
library("crch")
NGR_crch <- crch(temp ~ ensmean | enssd, data = temptrain)
library("gamlss")
NGR_gamlss <- gamlss(temp ~ ensmean, sigma.formula = ~ enssd, 
  data = temptrain)


###################################################
### code chunk number 25: NGR-fit2
###################################################
NGR_crch2 <- crch(temp ~ ensmean | I(enssd^2), data = temptrain, 
  link.scale = "quad", type = "crps")


###################################################
### code chunk number 26: NGR-fit3
###################################################
library("ensembleMOS")
NGR_ensembleMOS <- fitMOS(temptrain_eD, model = "normal")


###################################################
### code chunk number 27: NGR-coef
###################################################
coef_ensembleMOS <- c(NGR_ensembleMOS$a, 11 * NGR_ensembleMOS$B[1], 
  NGR_ensembleMOS$c, NGR_ensembleMOS$d)


###################################################
### code chunk number 28: NGR-coef2
###################################################
rbind(crch = coef(NGR_crch), 
  gamlss = c(coef(NGR_gamlss), coef(NGR_gamlss, what = "sigma")),
  ensembleMOS = coef_ensembleMOS, crch2 = coef(NGR_crch2))


###################################################
### code chunk number 29: BMA-fit
###################################################
BMA <- fitBMA(temptrain_eD, model = "normal")


###################################################
### code chunk number 30: ensdress-fit
###################################################
smuy2 <- var(MOS$residuals)
st2 <- mean(temptrain$enssd^2)
dress <- sqrt(smuy2 - (1 + 1/8) * st2)


###################################################
### code chunk number 31: ensdressvsBMA
###################################################
rbind(BMA = c(BMA$biasCoefs[,1], BMA$sd), 
  ensdress  = c(coef(MOS), sd = dress))


###################################################
### code chunk number 32: AKD-fit
###################################################
library("SpecsVerification")
(AKD <- FitAkdParameters(ens = as.matrix(temptrain[,2:12]), 
  obs = temptrain$temp))


###################################################
### code chunk number 33: NGRfc
###################################################
mean_NGR  <- predict(NGR_crch,  newdata = temptest, type = "location")
sd_NGR    <- predict(NGR_crch,  newdata = temptest, type = "scale")
mean_NGR2 <- predict(NGR_crch2, newdata = temptest, type = "location")
sd_NGR2   <- predict(NGR_crch2, newdata = temptest, type = "scale")


###################################################
### code chunk number 34: NGRfc2
###################################################
x <- seq(-10, 10, 0.1)
pdf_NGR   <- dnorm(x, mean_NGR[1], sd_NGR[1])
cdf_NGR   <- pnorm(0, mean_NGR, sd_NGR) 


###################################################
### code chunk number 35: NGRfc3
###################################################
quant_NGR <- predict(NGR_crch, newdata = temptest, type = "quantile", 
  at = c(0.25, 0.5, 0.75))


###################################################
### code chunk number 36: BMAfc
###################################################
corrected_BMA   <- apply(temptest[,2:12], 2, 
  function(x) BMA$biasCoefs[,1] %*% rbind(1, x))


###################################################
### code chunk number 37: ensdressfc
###################################################
corrected_dress <- apply(temptest[,2:12], 2, 
  function(x) coef(MOS) %*% rbind(1, x))


###################################################
### code chunk number 38: AKDfc
###################################################
AKDobj <- DressEnsemble(ens = as.matrix(temptest[2:12]), 
  dressing.method = "akd", parameters = as.list(AKD))


###################################################
### code chunk number 39: Ch11.Rnw:649-651
###################################################
dressobj <- list(ens = corrected_dress,
  ker.wd = matrix(dress, nrow = nrow(corrected_dress), ncol = 11))


###################################################
### code chunk number 40: BMA-ensdressfc
###################################################
cdf_BMA   <- cdf(BMA, temptest_eD, values = 0)
cdf_AKD   <- rowMeans(pnorm(0,   AKDobj$ens,   AKDobj$ker.wd))
cdf_dress <- rowMeans(pnorm(0, dressobj$ens, dressobj$ker.wd))


###################################################
### code chunk number 41: densityplot
###################################################
par(mfrow = c(1,4))
plot(x, pdf_NGR, type = "l",  xlab = "Temperature", ylab = "Density", 
  lwd = 3, main = "NGR")
abline(v = temptest$temp[1], col = "orange", lwd = 3)

plot(BMA, temptest_eD[1,])
title(main = "BMA")

plot(density(dressobj$ens[1,], bw = dressobj$ker.wd[1, 1]), 
  xlab= "Temperature",  main = "ensemble dressing", lwd = 3)
abline(v = temptest$temp[1], col = "orange", lwd = 3)


plot(density(AKDobj$ens[1,], bw = AKDobj$ker.wd[1, 1]), 
  xlab = "Temperature", main = "AKD", lwd = 3)
abline(v = temptest$temp[1], col = "orange", lwd = 3)


###################################################
### code chunk number 42: fig_densityplot
###################################################
par(mfrow = c(1,4))
plot(x, pdf_NGR, type = "l",  xlab = "Temperature", ylab = "Density", 
  lwd = 3, main = "NGR")
abline(v = temptest$temp[1], col = "orange", lwd = 3)

plot(BMA, temptest_eD[1,])
title(main = "BMA")

plot(density(dressobj$ens[1,], bw = dressobj$ker.wd[1, 1]), 
  xlab= "Temperature",  main = "ensemble dressing", lwd = 3)
abline(v = temptest$temp[1], col = "orange", lwd = 3)


plot(density(AKDobj$ens[1,], bw = AKDobj$ker.wd[1, 1]), 
  xlab = "Temperature", main = "AKD", lwd = 3)
abline(v = temptest$temp[1], col = "orange", lwd = 3)


###################################################
### code chunk number 43: quantilepred
###################################################
plot(quant_NGR[1:20, 2], type = "l", lty = 2, ylab = "2m temperature",
  xlab = "date", xaxt = "n", ylim = c(-15, 10))
axis(1, at = seq(1, 20, 6), temptest$date[seq(1, 20, 6)])
polygon(c(1:20, 20:1), c(quant_NGR[1:20, 1], quant_NGR[20:1, 3]), 
  col = gray(0.1, alpha = 0.1), border = FALSE)
lines(temptest$temp[1:20]) 


###################################################
### code chunk number 44: cdfpred
###################################################
plot(cdf_NGR[1:20], type = "l", ylab = "Pr(T<0)",
  xlab = "date", xaxt = "n", ylim = c(0, 1))
axis(1, at = seq(1, 20, 6), temptest$date[seq(1, 20, 6)])
points(temptest$temp<0) 


###################################################
### code chunk number 45: fig_quantilepred
###################################################
par(mfrow = c(1,2))
plot(quant_NGR[1:20, 2], type = "l", lty = 2, ylab = "2m temperature",
  xlab = "date", xaxt = "n", ylim = c(-15, 10))
axis(1, at = seq(1, 20, 6), temptest$date[seq(1, 20, 6)])
polygon(c(1:20, 20:1), c(quant_NGR[1:20, 1], quant_NGR[20:1, 3]), 
  col = gray(0.1, alpha = 0.1), border = FALSE)
lines(temptest$temp[1:20]) 
plot(cdf_NGR[1:20], type = "l", ylab = "Pr(T<0)",
  xlab = "date", xaxt = "n", ylim = c(0, 1))
axis(1, at = seq(1, 20, 6), temptest$date[seq(1, 20, 6)])
points(temptest$temp<0) 


###################################################
### code chunk number 46: crps
###################################################
library("scoringRules")
crps_all <- cbind(
  NGR1  = scoringRules::crps(temptest$temp, family = "normal", 
    mean = mean_NGR, sd = sd_NGR),
  NGR2  = scoringRules::crps(temptest$temp, family = "normal", 
    mean = mean_NGR2, sd = sd_NGR2),
  BMA   = ensembleBMA::crps(BMA, temptest_eD)[, 2],
  dress = DressCrps(dressobj, temptest$temp),
  AKD   = DressCrps(AKDobj, temptest$temp))


###################################################
### code chunk number 47: bootmean
###################################################
bootmean <- function(scores, nsamples = 250) {
  boot <- NULL
  for(i in 1:nsamples) { 
    bindex <- sample(nrow(scores), replace = TRUE)
    boot <- rbind(boot, colMeans(scores[bindex,]))
  }
  boot
}


###################################################
### code chunk number 48: crps2
###################################################
boxplot(bootmean(crps_all), ylab = "CRPS")


###################################################
### code chunk number 49: ign
###################################################
ign_all <- cbind(
  NGR1  = -dnorm(temptest$temp, mean_NGR, sd_NGR, log = TRUE),
  NGR2  = -dnorm(temptest$temp, mean_NGR2, sd_NGR2, log = TRUE),
  BMA   = -log(rowSums(BMA$weights * dnorm(temptest$temp, corrected_BMA, BMA$sd))),
  dress = -log(rowMeans(dnorm(temptest$temp, dressobj$ens, dressobj$ker.wd))),
  AKD   = -log(rowMeans(dnorm(temptest$temp, AKDobj$ens, AKDobj$ker.wd))))


###################################################
### code chunk number 50: ign2
###################################################
boxplot(bootmean(ign_all), ylab = "IS")


###################################################
### code chunk number 51: fig_bootscores
###################################################
par(mfrow=c(1,2), mar=c(4,4,2,2))
boxplot(bootmean(crps_all), ylab = "CRPS")
boxplot(bootmean(ign_all), ylab = "IS")


###################################################
### code chunk number 52: pit
###################################################
pit <- cbind(
  NGR1  = pnorm(temptest$temp, mean_NGR,  sd_NGR),
  NGR2  = pnorm(temptest$temp, mean_NGR2, sd_NGR2),
  BMA   = pit(BMA, temptest_eD),
  dress = rowMeans(pnorm(temptest$temp, dressobj$ens, dressobj$ker.wd)),
  AKD   = rowMeans(pnorm(temptest$temp, AKDobj$ens, AKDobj$ker.wd)))


###################################################
### code chunk number 53: pithist
###################################################
par(mfrow = c(1, ncol(pit)))
for(model in colnames(pit)){
  hist(pit[, model], main = model, freq = FALSE, 
    xlab = "", ylab = "", ylim = c(0, 1.6))
  abline(h = 1, lty = 2)
}


###################################################
### code chunk number 54: fig_pithist
###################################################
par(mar=c(2,2,2,2))
par(mfrow = c(1, ncol(pit)))
for(model in colnames(pit)){
  hist(pit[, model], main = model, freq = FALSE, 
    xlab = "", ylab = "", ylim = c(0, 1.6))
  abline(h = 1, lty = 2)
}


###################################################
### code chunk number 55: data-rain
###################################################
data("rain", package = "ensemblepp")


###################################################
### code chunk number 56: data-rain2
###################################################
dim(rain)
names(rain)


###################################################
### code chunk number 57: data-rain3
###################################################
rain <- sqrt(rain)


###################################################
### code chunk number 58: data-rain4
###################################################
rain$date <- as.Date(rownames(rain))
rain <- rain[format(rain$date, "%m") %in% c("12", "01", "02"),]


###################################################
### code chunk number 59: data-rain5
###################################################
rain$ensmean <- apply(rain[,2:12], 1, mean)
rain$enssd   <- apply(rain[,2:12], 1, sd)


raintrain <- rain[rain$date < "2010-03-01",]
raintest  <- rain[rain$date > "2010-03-01",]


###################################################
### code chunk number 60: rainplots
###################################################
par(mfrow = c(2,2))
plot(rain~ensmean, raintrain, col = gray(0.2, alpha = 0.4),
  main = "Scatterplot")
abline(0, 1, lty = 2)

rank <- apply(raintrain[,1:12], 1, rank)[1,]
hist(rank, breaks = 0:12 + 0.5, main = "Verification Rank Histogram")

sdcat <- cut(raintrain$enssd, quantile(raintrain$enssd, seq(0, 1, 0.2)))
boxplot(abs(rain-ensmean)~sdcat, raintrain, ylab = "absolute error",
  xlab = "ensemble standard deviation", main = "Spread-Skill")

hist(rain$rain, xlab = "square root of precipitation", main = "Histogram")


###################################################
### code chunk number 61: fig_rainplots
###################################################
par(mfrow = c(2,2))
plot(rain~ensmean, raintrain, col = gray(0.2, alpha = 0.4),
  main = "Scatterplot")
abline(0, 1, lty = 2)

rank <- apply(raintrain[,1:12], 1, rank)[1,]
hist(rank, breaks = 0:12 + 0.5, main = "Verification Rank Histogram")

sdcat <- cut(raintrain$enssd, quantile(raintrain$enssd, seq(0, 1, 0.2)))
boxplot(abs(rain-ensmean)~sdcat, raintrain, ylab = "absolute error",
  xlab = "ensemble standard deviation", main = "Spread-Skill")

hist(rain$rain, xlab = "square root of precipitation", main = "Histogram")


###################################################
### code chunk number 62: NGRrainfit
###################################################
cNLR <- crch(rain ~ ensmean | enssd, data = raintrain, left = 0, 
  dist = "logistic")


###################################################
### code chunk number 63: data-rain6
###################################################
raintrain_eD <- ensembleData(forecasts = raintrain[,2:12],
  dates = raintrain$date, observations = raintrain$rain,
  forecastHour = 24, initializationTime = "00", exchangeable = rep(1, 11))
raintest_eD <- ensembleData(forecasts = raintest[,2:12],
  dates = raintest$date, observations = raintest$rain,
  forecastHour = 24, initializationTime = "00", exchangeable = rep(1, 11))


###################################################
### code chunk number 64: BMArainfit
###################################################
gBMA <- fitBMA(raintrain_eD, model = "gamma0", 
  control = controlBMAgamma0(power = 1))


###################################################
### code chunk number 65: logregrainfit
###################################################
logreg <- glm(rain > 0 ~ ensmean, data = raintrain, family = binomial())


###################################################
### code chunk number 66: hlogregrainfit
###################################################
library("glmx")
hlogreg <- hetglm(rain > 0 ~ ensmean | enssd, data = raintrain, 
  family = binomial())


###################################################
### code chunk number 67: rain-quantiles
###################################################
q <- unique(quantile(raintrain$rain, seq(0.1, 0.9, 0.1)))


###################################################
### code chunk number 68: combinelogreg
###################################################
logreg2 <- hetlogreg2 <- list()
for(i in 1:length(q)){
  logreg2[[i]] <- glm(rain <= q[i] ~ ensmean, data = raintrain, 
    family = binomial())
  hetlogreg2[[i]] <- hetglm(rain <= q[i] ~ ensmean | enssd, 
    data = raintrain, family = binomial())
}


###################################################
### code chunk number 69: HXLRfit
###################################################
HXLR <- hxlr(rain ~ ensmean | enssd, data = raintrain, thresholds = q)


###################################################
### code chunk number 70: raincat
###################################################
raintrain$raincat <- cut(raintrain$rain, c(-Inf, q, Inf)) 


###################################################
### code chunk number 71: clmrainfit
###################################################
library("ordinal")
OLR <- clm(raincat ~ ensmean, scale = ~ enssd, data = raintrain)


###################################################
### code chunk number 72: NGRrainfc
###################################################
location_cNLR <- predict(cNLR, newdata = raintest, type = "location")
scale_cNLR    <- predict(cNLR, newdata = raintest, type = "scale")


###################################################
### code chunk number 73: quantilepred2
###################################################
par(mfrow = c(1,2))
quant_cNLR <- predict(cNLR, newdata = raintest, type = "quantile", 
  at = c(0.25, 0.5, 0.75))
plot(quant_cNLR[1:20, 2], type = "l", lty = 2, ylim = c(0, 5),
  ylab = "square root of precipitation", xlab = "date", xaxt = "n")
axis(1, at = seq(1, 20, 6), raintest$date[seq(1, 20, 6)])
polygon(c(1:20, 20:1), c(quant_cNLR[1:20, 1], quant_cNLR[20:1, 3]), 
  col = gray(0.1, alpha = 0.1), border = FALSE)
lines(raintest$rain[1:20]) 


###################################################
### code chunk number 74: cdfpred2
###################################################
cdf_cNLR <- sapply(q, function(q) plogis(q, location_cNLR, scale_cNLR))
plot(1 - cdf_cNLR[1:20], type = "l", ylab = "Pr(rain>0)",
  xlab = "date", xaxt = "n", ylim = c(0,1))
axis(1, at = seq(1, 20, 6), temptest$date[seq(1, 20, 6)])
points(raintest$rain>0) 


###################################################
### code chunk number 75: fig_quantilepred2
###################################################
par(mfrow = c(1,2))
quant_cNLR <- predict(cNLR, newdata = raintest, type = "quantile", 
  at = c(0.25, 0.5, 0.75))
plot(quant_cNLR[1:20, 2], type = "l", lty = 2, ylim = c(0, 5),
  ylab = "square root of precipitation", xlab = "date", xaxt = "n")
axis(1, at = seq(1, 20, 6), raintest$date[seq(1, 20, 6)])
polygon(c(1:20, 20:1), c(quant_cNLR[1:20, 1], quant_cNLR[20:1, 3]), 
  col = gray(0.1, alpha = 0.1), border = FALSE)
lines(raintest$rain[1:20]) 
cdf_cNLR <- sapply(q, function(q) plogis(q, location_cNLR, scale_cNLR))
plot(1 - cdf_cNLR[1:20], type = "l", ylab = "Pr(rain>0)",
  xlab = "date", xaxt = "n", ylim = c(0,1))
axis(1, at = seq(1, 20, 6), temptest$date[seq(1, 20, 6)])
points(raintest$rain>0) 


###################################################
### code chunk number 76: raindensities
###################################################
par(mfrow = c(1,3))
x <- c(0, seq(1e-8, 6.5, 0.1))
## cNLR
plot(x, dclogis(x, location_cNLR[2], scale_cNLR[2], left = 0), type = "l",
  lwd = 3, xlab = "Square root of precipitation", ylab = "PDF", main = "cNLR")
abline(v = raintest$rain[2], lwd = 3, col = "orange")
## BMA
plot(gBMA, raintest_eD[2,])
title(main = "BMA")
## HXLR
location_HXLR <- predict(HXLR, newdata = raintest, type = "location")
scale_HXLR    <- predict(HXLR, newdata = raintest, type = "scale")
plot(x, dclogis(x, location_HXLR[2], scale_HXLR[2], left = 0), type = "l",
  lwd = 3, xlab = "Square root of precipitation", ylab = "PDF", main = "HXLR")
abline(v = raintest$rain[2], lwd = 3, col = "orange")


###################################################
### code chunk number 77: fig_raindensities
###################################################
par(mfrow = c(1,3))
x <- c(0, seq(1e-8, 6.5, 0.1))
## cNLR
plot(x, dclogis(x, location_cNLR[2], scale_cNLR[2], left = 0), type = "l",
  lwd = 3, xlab = "Square root of precipitation", ylab = "PDF", main = "cNLR")
abline(v = raintest$rain[2], lwd = 3, col = "orange")
## BMA
plot(gBMA, raintest_eD[2,])
title(main = "BMA")
## HXLR
location_HXLR <- predict(HXLR, newdata = raintest, type = "location")
scale_HXLR    <- predict(HXLR, newdata = raintest, type = "scale")
plot(x, dclogis(x, location_HXLR[2], scale_HXLR[2], left = 0), type = "l",
  lwd = 3, xlab = "Square root of precipitation", ylab = "PDF", main = "HXLR")
abline(v = raintest$rain[2], lwd = 3, col = "orange")


###################################################
### code chunk number 78: BMArainfc
###################################################
cdf_gBMA <- cdf(gBMA, raintest_eD, values = q)


###################################################
### code chunk number 79: logregrainfc
###################################################
cdf_logreg <- sapply(logreg2, function(mod)
  predict(mod, newdata = raintest, type = "response"))
cdf_hetlogreg <- sapply(hetlogreg2, function(mod) 
  predict(mod, newdata = raintest, type = "response"))


###################################################
### code chunk number 80: HXLR-OLRrainfc
###################################################
cdf_HXLR <- predict(HXLR, newdata = raintest, type = "cumprob")
cdf_OLR <- predict(OLR, newdata = raintest, type = "cum.prob")$cprob1[,-7]


###################################################
### code chunk number 81: logregrainfc
###################################################
CDF <- list(cNLR = cdf_cNLR, BMA = cdf_gBMA, logreg = cdf_logreg, 
  hlogreg = cdf_hetlogreg, HXLR = cdf_HXLR, OLR = cdf_OLR)


###################################################
### code chunk number 82: brierrain
###################################################
brier_all <- NULL
for(n in names(CDF)) {
  brier_all <- cbind(brier_all, ((raintest$rain <= 0) - CDF[[n]][,1])^2)
}
colnames(brier_all) <- names(CDF)


###################################################
### code chunk number 83: brier
###################################################
boxplot(bootmean(brier_all), las = 2, ylab = "Brier score")


###################################################
### code chunk number 84: rpsprep
###################################################
cdf_obs <- sapply(q, function(q) raintest$rain <= q)
rps_all <- NULL
for(n in names(CDF)) {
  rps_all <- cbind(rps_all, 
    rowSums((CDF[[n]] - cdf_obs)^2)/(ncol(cdf_HXLR) - 1))
}
colnames(rps_all) <- names(CDF)


###################################################
### code chunk number 85: rps
###################################################
boxplot(bootmean(rps_all), las = 2, ylab = "RPS")


###################################################
### code chunk number 86: fig_bootscores2
###################################################
par(mfrow=c(1,2), mar=c(4.5,4,2,2))
boxplot(bootmean(brier_all), las = 2, ylab = "Brier score")
boxplot(bootmean(rps_all), las = 2, ylab = "RPS")


###################################################
### code chunk number 87: brierdecomp
###################################################
sapply(CDF, function(x) BrierDecomp(x[,1], y = (raintest$rain <= 0))[1,])


###################################################
### code chunk number 88: reliability
###################################################
par(mfrow=c(2,3))
for(n in 1:length(names(CDF))){
  ReliabilityDiagram(1 - CDF[[n]][,1], (raintest$rain > 0), plot=TRUE)
  par(mfg = c((n-1) %/% 3 + 1, (n-1) %% 3 + 1))
  title(main = names(CDF)[n])
}


###################################################
### code chunk number 89: fig_reliability
###################################################
par(mfrow=c(2,3))
for(n in 1:length(names(CDF))){
  ReliabilityDiagram(1 - CDF[[n]][,1], (raintest$rain > 0), plot=TRUE)
  par(mfg = c((n-1) %/% 3 + 1, (n-1) %% 3 + 1))
  title(main = names(CDF)[n])
}


###################################################
### code chunk number 90: roc
###################################################
library("pROC")
par(mfrow = c(2, 3))
for(n in names(CDF)){
  rocplot <- roc((raintest$rain > 0) ~ I(1 - CDF[[n]][, 1]), plot=TRUE,
    main = n)
  text(0.2, 0.2, paste("AUC =", round(rocplot$auc, digits = 4)))
}


###################################################
### code chunk number 91: fig_roc
###################################################
library("pROC")
par(mfrow = c(2, 3))
for(n in names(CDF)){
  rocplot <- roc((raintest$rain > 0) ~ I(1 - CDF[[n]][, 1]), plot=TRUE,
    main = n)
  text(0.2, 0.2, paste("AUC =", round(rocplot$auc, digits = 4)))
}


###################################################
### code chunk number 92: rpsprep2
###################################################
cdf_obs <- sapply(q, function(q) raintest$rain <= q)
rps_all <- NULL
for(n in names(CDF)) {
  rps_all <- cbind(rps_all, 
    rowSums((CDF[[n]] - cdf_obs)^2)/(ncol(cdf_HXLR) - 1))
}
colnames(rps_all) <- names(CDF)


###################################################
### code chunk number 93: rps2
###################################################
boxplot(bootmean(rps_all), las = 2, ylab = "RPS")


###################################################
### code chunk number 94: temprain
###################################################
plot(raintrain$rain~temptrain$temp, xlab = "temperature", 
  ylab = "precipitation")
abline(lm(raintrain$rain~temptrain$temp))


###################################################
### code chunk number 95: fig_raintemp
###################################################
plot(raintrain$rain~temptrain$temp, xlab = "temperature", 
  ylab = "precipitation")
abline(lm(raintrain$rain~temptrain$temp))


###################################################
### code chunk number 96: tr-temp-rain
###################################################
trtemp <- qnorm(pnorm(temptrain$temp, 
  fitted(NGR_crch, type = "location"), 
  fitted(NGR_crch, type = "scale"))) 
trrain <- qnorm(clogispit(raintrain$rain, 
  fitted(cNLR, type = "location"), 
  fitted(cNLR, type = "scale")))


###################################################
### code chunk number 97: cov-temp-rain
###################################################
covmatrix <- cov(cbind(trtemp, trrain))


###################################################
### code chunk number 98: nscen
###################################################
nscen <- 11


###################################################
### code chunk number 99: MVfc
###################################################
library("mvtnorm")
GCtemp <- GCrain <- NULL
for(i in 1:nrow(temptest)) {
  sim <- rmvnorm(nscen, c(0, 0), sigma = covmatrix)
  GCtemp <- rbind(GCtemp, qnorm(pnorm(sim[,1]), mean_NGR[i], sd_NGR[i]))
  GCrain <- rbind(GCrain, 
    qclogis(pnorm(sim[,2]), location_cNLR[i], scale_cNLR[i], left = 0))
}


###################################################
### code chunk number 100: MVfc2
###################################################
UNIVtemp <- UNIVrain <- NULL
for(i in 1:nrow(temptest)) {
  UNIVtemp <- rbind(UNIVtemp, rnorm(nscen, mean_NGR[i], sd_NGR[i]))
  UNIVrain <- rbind(UNIVrain, 
    rclogis(nscen, location_cNLR[i], scale_cNLR[i], left = 0))
}


###################################################
### code chunk number 101: MVfc3
###################################################
ECCtemp <- ECCrain <- NULL
for(i in 1:nrow(temptest)) {
  ECCtemp <- rbind(ECCtemp, sort(UNIVtemp[i,])[rank(temptest[i, 2:12])])
  ECCrain <- rbind(ECCrain, sort(UNIVrain[i,])[rank(raintest[i, 2:12])])
}


###################################################
### code chunk number 102: MVfc4
###################################################
SStemp <- SSrain <- NULL
for(i in 1:nrow(temptest)) {
  ind <- sample(nrow(temptrain), nscen)
  SStemp <- rbind(SStemp, sort(UNIVtemp[i,])[rank(temptrain$temp[ind])])
  SSrain <- rbind(SSrain, sort(UNIVrain[i,])[rank(raintrain$rain[ind])])
}


###################################################
### code chunk number 103: MVverification
###################################################
es <- vs <- NULL
for(i in 1:nrow(temptest)) {
  obs <- c(temptest$temp[i], raintest$rain[i])  
  es <- rbind(es, c(
    UNIV = es_sample(obs, rbind(UNIVtemp[i,], UNIVrain[i,])),
    ECC  = es_sample(obs, rbind(ECCtemp[i,], ECCrain[i,])),
    GC   = es_sample(obs, rbind(GCtemp[i,], GCrain[i,])),
    SS   = es_sample(obs, rbind(SStemp[i,], SSrain[i,]))))
  vs <- rbind(vs, c(
    UNIV = vs_sample(obs, rbind(UNIVtemp[i,], UNIVrain[i,])),
    ECC  = vs_sample(obs, rbind(ECCtemp[i,], ECCrain[i,])),
    GC   = vs_sample(obs, rbind(GCtemp[i,], GCrain[i,])),
    SS   = vs_sample(obs, rbind(SStemp[i,], SSrain[i,]))))

}


###################################################
### code chunk number 104: bootscores3
###################################################
par(mfrow = c(1, 2))
boxplot(bootmean(es), ylab = "ES")
boxplot(bootmean(vs), ylab = "VS")


###################################################
### code chunk number 105: fig_bootscores3
###################################################
par(mfrow = c(1, 2))
boxplot(bootmean(es), ylab = "ES")
boxplot(bootmean(vs), ylab = "VS")


###################################################
### code chunk number 106: mvrank
###################################################
r <- NULL
for(i in 1:nrow(temptest)) {
  obs <- c(temptest$temp[i], raintest$rain[i])
  r <- rbind(r, c(
    UNIV = mv.rank(cbind(obs, rbind(UNIVtemp[i,], UNIVrain[i,])))[1],
    ECC  = mv.rank(cbind(obs, rbind(ECCtemp[i,], ECCrain[i,])))[1],
    GC   = mv.rank(cbind(obs, rbind(GCtemp[i,], GCrain[i,])))[1],
    SS   = mv.rank(cbind(obs, rbind(SStemp[i,], SSrain[i,])))[1]))
}


###################################################
### code chunk number 107: mvrankhist
###################################################
par(mfrow = c(1, 4))
for(n in colnames(r)){
  hist(r[,n], freq = FALSE, breaks = 0:12 + 0.5, main = n, ylim = c(0, 0.12))
  abline(h = 1/11, lty = 2)
}


###################################################
### code chunk number 108: fig_mvrankhist
###################################################
par(mfrow = c(1, 4))
for(n in colnames(r)){
  hist(r[,n], freq = FALSE, breaks = 0:12 + 0.5, main = n, ylim = c(0, 0.12))
  abline(h = 1/11, lty = 2)
}


###################################################
### code chunk number 109: Ch11.Rnw:1740-1741 (eval = FALSE)
###################################################
## clogispit <- function(obs, location, scale, left = 0) {
##   pit <- pclogis(obs, location, scale, left = left)
##   pit[obs <= 0] <- runif(sum(obs <= 0), 0, pit[obs <= 0])
##   return(pit)
## }


###################################################
### code chunk number 110: app:mvrank (eval = FALSE)
###################################################
## mv.rank <- function(x) {
##   d <- dim(x)
##   x.prerank <- numeric(d[2])
##   for(i in 1:d[2]) {
##     x.prerank[i] <- sum(apply(x<=x[,i],2,all))
##   }
##   x.rank <- rank(x.prerank,ties="random")
##   return(x.rank)
## }


