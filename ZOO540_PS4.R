# Angelica, Biz, Quinn
# PS4
# due 17 October, 2023 by midnight
library(lme4)
library(lmerTest)
# The mvtnorm package has a function for randomly simulating from a multivariate normal distribution.
library(mvtnorm)
# This problem set is based on correlated_data Ch2 ("Ch_2_10Aug18.R").

# For the homework, ONLY TURN IN THIS FILE WITH THE CODE THAT YOU ADDED TO IT. Start with the code from correlated_data ("Ch_1_10Aug18.R") that is below and add to it anything you need. Identify you new code (so I can find it) by placing it between marker rows #~~~~~~~~~~~~~~~~~~~~~~~~~~~. 

# 1. What questions do you have about the material in chapter 2? What needs more explanation? I'm serious about asking this question, because I want to improve the book. (NOTE: This requires NO NEW R CODE.)

# 2. As discussed in subsection 2.6.5, the LM should give good type I errors even when applied to binary data, provided the residuals are homoscedastic. This will be true for the simplest case of only a single predictor (independent) variable, because under the null hypothesis H0:b1 = 0, the residuals will be homoscedastic (since the predicted values for all points are the same). However, if there is a second predictor variable that is correlated to the first, then under the null hypothesis H0:b1 = 0 the residuals will not be homoscedastic. This could generate incorrect type I errors for the LM. In the simulations of binary data, this didn't seem to cause problems with type I error for the LM (section 2.6). For this exercise, investigate this problem for binomial data with more than two (i.e., 0/1) outcomes by computing the rate at which the null hypothesis H0:b1=0 is rejected in simulated data under the null hypothesis, exploring the case when there is a second independent variable x2 that is highly correlated with x1. You can modify the code from section 2.6 for this.
#~~~~~~~~~~~~~~~~~~~~~~~~~~
n <- 1000
b0 <- 2
b1 <- 0
b2 <- 1.5
var.x <- 1
cov.x <- .9

size <- 2

x <- rmvnorm(n, sigma = matrix(c(var.x, cov.x, cov.x, var.x), nrow=2))
x1 <- x[,1]
x2 <- x[,2]
Y <- rbinom(n=n, size=size, prob=inv.logit(b0 + b1*x1 + b2*x2))
dat <- data.frame(x1=x1, x2=x2, Y=Y)

# Fit the reduced LM model
mod.reduced <- lm(Y ~ x2, data=dat)
dat$r <- mod.reduced$resid

# Fig. 2.7
# Plot the residuals
par(mfrow=c(1,3))
plot(r ~ x1, data=dat, xlab="Values of x1", ylab="Residuals")

# Bin and plot the variances of the residuals
dat <- dat[order(dat$x1),]
bins <- round(n*(0:10)/10)
resid <- data.frame(x1=bins[-1], var=NA)
for(i in 2:11){
  resid$x1[i-1] <- mean(dat$x1[bins[i-1]:bins[i]])
  resid$var[i-1] <- var(dat$r[bins[i-1]:bins[i]])
}
plot(var ~ x1, data=resid, ylim=c(0,max(var)), xlab="Binned values of x1", ylab="Residual variance")

plot(Y~x1, data=dat)


n.list <- c(10,13,20,30,50,100,200)
nsims <- 2000

reject <- data.frame(n=n.list, lm=NA)
i.n <- 0 
for(n in n.list){
  i.n <- i.n +1
  Pvalues <- data.frame(glm.Wald=array(NA, dim=nsims), glm.con=NA, glm.LRT=NA, glm.boot0=NA, lm=NA)
  for(i in 1:nsims){
  x <- rmvnorm(n, sigma = matrix(c(var.x, cov.x, cov.x, var.x), nrow=2))
  x1 <- x[,1]
  x2 <- x[,2]
  Z <- b0+b1*x1 + b2*x2
  prob <- inv.logit(Z)
  Y <- rbinom(n=n, size=1, prob=inv.logit(b0 + b1*x1 + b2*x2))
  mod.lm <- lm(Y ~ x1+x2)
  Pvalues$lm[i] <- summary(mod.lm)$coef[2,4]
  
  reject$lm[i.n] <- mean(Pvalues$lm , 0.05, na.rm = T)
}
}
round(reject, digits=3)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. I illustrated the two bootstraps (bootstrapping the estimator of b1 and bootstrapping under the null hypothesis) using a GLM simulation and model. Analyze these bootstraps for a LM. Specifically, create a version of Fig. 2.6 for a LM in which you simulate data using a LM, and you then fit the simulated data with an LM. Both bootstraps give you a distribution of bootstrap estimates: for the bootstrap of b1, the estimates approximate the estimator of b1, and for the bootstrap under the null hypothesis, the estimates approximate the distribution of the deviance. For a LM, the estimator of b1 should be t-distributed, and the estimator of the deviance should be chi-squared distributed. By plotting the t-distribution with the bootstrapped estimator of b1 (like the blue line in the left panel of Fig. 2.6) and the chi-squared distribution with the bootstrapped estimator of the deviance (like the blue line in the right panel of Fig. 2.6) , confirm that this is true.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fit full and reduced models to the data.
mod.dat <- lm(Y ~ x1 + x2, data=dat)
mod.reduced <- lm(Y ~ x2, data=dat)
dev.dat <- 2*(logLik(mod.dat) - logLik(mod.reduced))[1]

# Simulate the data using the reduced model.
boot0 <- data.frame(LLR=rep(NA, nboot), converge=NA)
dat.boot <- dat
for(i in 1:nboot){
  dat.boot$Y <- simulate(mod.reduced)[[1]]
  mod.lm <- update(mod.dat, data=dat.boot)
  mod.lm0 <- update(mod.reduced, data=dat.boot)
  boot0$dev[i] <- 2*(logLik(mod.lm) - logLik(mod.lm0))[1]
#  boot0$converge[i] <- mod.lm$converge
}
# Remove the cases when glm() did not converge
#boot0 <- boot0[boot0$converge == T,]
#pvalue <- mean(boot0$dev > dev.dat)
#pvalue

# This is the p-value from a LRT
pchisq(dev.dat, df=1, lower.tail=F)

# Fig. 2.6 
# Comparing both bootstraps
par(mfrow=c(1,2))

# For small n, it might be necessary to remove poorly fit simulations
boot <- boot[boot$b1 < (mean(boot$b1) + 2.5*sd(boot$b1)),]

# Plot of parametric bootstrap
hist(boot$b1, xlab="Bootstrap of b1", main = paste("b1 = ", round(mod.dat$coef[2],2), ",  Estimate mean = ", round(mean(boot$b1),2), sep=''), freq=F, breaks=40)
#lines(c(mod.dat$coef[2],mod.dat$coef[2]), c(0,2), col="red")
#lines(c(0,0), c(0,2), col="green")
#lines(.1*(-100:100), dnorm(.1*(-100:100), mean=mod.dat$coef[2], sd=summary(mod.dat)$coef[2,2]), col="blue")

# Plot of parametric bootstrap around H0
hist(boot0$dev, xlab="Bootstrap of deviance", main = paste("Deviance = ", round(dev.dat,2), sep=''), freq=F, breaks=40)
#lines(c(dev.dat,dev.dat), c(0,1), col="red")
lines(.1*(1:100), dchisq(.1*(1:100), df=1), col="blue")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. For the bootstrap estimator of b1, the P-value is given by

# pvalue <- 2*min(mean(boot$b1 < 0), mean(boot$b1 > 0))

# while for the bootstrap estimator of the deviance, the P-value is given by 

# pvalue <- mean(boot0$LLR > LLR.dat)

# Explain why these equations are different. Specifically, why does the P-value for the bootstrap of the estimator of b1 take the minimum of two values? (You DON'T need to write any R code.)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#The two equations are different because they address different hypotheses and test different aspects of the model. The P-value for the estimator of b1 considers both sides of the distribution (two-tailed test), so the minimum operator is used to account for both cases, while the P-value for the deviance estimator is based on a one-tailed test, so there's no need for the minimum operator.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. When analyzing the binary regression model, I left out one of the best methods: logistic regression using a Firth correction. I won't go into details to explain how logistic regression with a Firth correction works, but the basic idea is to fit the regression using ML while penalizing the likelihood so that it doesn't show (as much) bias as the binomial glm. To see how logistic regression with a Firth correction performs, use the function logistf() in the package {logistf} to produce the same power figure as in Fig. 2.8 (subsection 2.6.6). Don't bother including the bootstrap tests, since they take more time: just use the Wald and LRT tests with the GLM fits, and the standard (t-distribution) test for the LM. You can still compare how well logistf() does compared to the other methods, since the LM performed the same as the bootstrap of H0. For the simulations, you fit logistf() using the syntax

# mod.logistf <- logistf(Y ~ x1 + x2)

# You can then extract the P-value for x1 with mod.logistf$prob[2].
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(logistf)

# Compute the rejection rates for the methods except the parametric bootstrap of b1 (which is done in the following code). Note that you will need to calculate the distribution of LLR for bootstrap of H0 using the code at the top of section "2.6 P-values for binary data" above.
n <- 30
b0 <- 1
b1 <- 0
b2 <- 1.5
var.x <- 1
cov.x <- 0
b1.list <- c(0, .2, .4, .6, .8, 1)
nsims <- 1000

reject <- data.frame(b1=b1.list, glm.Wald=NA, glm.Wald.con=NA, glm.LRT=NA, glm.boot0=NA, glm.boot0.con=NA, lm=NA, glm.con=NA, logistf=NA)
i.b1 <- 0
for(b1 in b1.list){
  i.b1 <- i.b1 + 1
  Pvalues <- data.frame(glm.Wald=array(NA,dim=nsims), glm.con=NA, glm.LRT=NA, glm.boot0=NA, lm=NA, logistf=NA)
  for(i in 1:nsims){
    x <- rmvnorm(n, sigma = matrix(c(var.x, cov.x, cov.x, var.x), nrow=2))
    x1 <- x[,1]
    x2 <- x[,2]
    Y <- rbinom(n=n, size=1, prob=inv.logit(b0 + b1*x1 + b2*x2))
    
    mod.glm <- glm(Y ~ x1 + x2, family = binomial)
    mod.glm0 <- glm(Y ~ x2, family = binomial)
    mod.lm <- lm(Y ~ x1 + x2)
    mod.logistf <- logistf(Y ~ x1 + x2)
    
    Pvalues$glm.Wald[i] <- summary(mod.glm)$coef[2,4]
    Pvalues$glm.LRT[i] <- pchisq(2*(logLik(mod.glm) - logLik(mod.glm0)), df=1, lower.tail=F)
    Pvalues$glm.deviance[i] <- 2*(logLik(mod.glm) - logLik(mod.glm0))
    Pvalues$lm[i] <- summary(mod.lm)$coef[2,4]
    
    Pvalues$glm.con[i] <- mod.glm$converge
    Pvalues$logistf[i] <-mod.logistf$prob[2]
  }
  reject$glm.Wald[i.b1] <- mean(Pvalues$glm.Wald < 0.05)
  reject$lm[i.b1] <- mean(Pvalues$lm < 0.05)
  reject$glm.Wald.con[i.b1] <- mean(Pvalues$glm.Wald[Pvalues$glm.con==T] < 0.05)
  reject$glm.LRT[i.b1] <- mean(Pvalues$glm.LRT < 0.05)
  reject$glm.boot0[i.b1] <- mean(Pvalues$glm.deviance > LLR.crit[LLR.crit$n == n,3])
  reject$glm.boot0.con[i.b1] <- mean(Pvalues$glm.deviance[Pvalues$glm.con==T] > LLR.crit[LLR.crit$n == n,3])
  reject$glm.con[i.b1] <- mean(Pvalues$glm.con)
  reject$logistf[i.b1] <- mean(Pvalues$logistf < 0.05)
  show(reject[i.b1,])
}
round(reject, digits=3)
# b1 glm.Wald glm.Wald.con glm.LRT glm.boot0 glm.boot0.con     lm glm.con
# 1 0.0   0.0228   0.02290997  0.0755    0.0536    0.05144695 0.0505  0.9952
# 2 0.2   0.0389   0.03912694  0.1001    0.0734    0.07071012 0.0723  0.9942
# 3 0.4   0.0715   0.07190989  0.1620    0.1238    0.12038620 0.1244  0.9943
# 4 0.6   0.1406   0.14160540  0.2748    0.2236    0.21945815 0.2234  0.9929
# 5 0.8   0.2431   0.24478904  0.4064    0.3407    0.33732756 0.3427  0.9931
# 6 1.0   0.3526   0.35576632  0.5500    0.4819    0.47845828 0.4820  0.9911

# Parametric bootstrap of b1. Note that this takes a lot of time to run.
#nsims <- 50
#nboot <- 20

#reject.boot <- data.frame(b1=b1.list, glm.boot=NA)
#i.b1 <- 0
#for(b1 in b1.list){
#  i.b1 <- i.b1 + 1
#  Pvalues <- data.frame(glm.boot=rep(NA, nsims))
#  for(i in 1:nsims){
#    x <- rmvnorm(n, sigma = matrix(c(var.x, cov.x, cov.x, var.x), nrow=2))
#    x1 <- x[,1]
#    x2 <- x[,2]
#    Y <- rbinom(n=n, size=1, prob=inv.logit(b0 + b1*x1 + b2*x2))
#    dat <- data.frame(x1=x1, x2=x2, Y=Y)
#    
#    mod.real <- glm(Y ~ x1 + x2, family = binomial, dat)
    # Parametric bootstrap
#    estimate <- data.frame(b1=array(NA,dim=nboot))
#    dat.boot <- dat
#    for(j in 1:nboot){
#      dat.boot$Y <- simulate(mod.real)[[1]]
#      mod.boot <- glm(Y ~ x1 + x2, family = binomial, data=dat.boot)
#      if(mod.boot$converge == T) estimate$b1[j] <- mod.boot$coef[2]
#    }
#    Pvalues$glm.boot[i] <- 2*min(mean(estimate$b1 < 0, na.rm=T), mean(estimate$b1 > 0, na.rm=T))
#  }
#  reject.boot$glm.boot[i.b1] <- mean(Pvalues$glm.boot < 0.05, na.rm=T)
#  show(reject.boot[i.b1,])
#}
# b1  glm.boot
# 1 0.0 0.1022044
# 2 0.2 0.1064257
# 3 0.4 0.1943888
# 4 0.6 0.2897384
# 5 0.8 0.4567404
# 6 1.0 0.5688259

#reject$glm.boot <- reject.boot$glm.boot
# Because the code takes so long to run, I saved it.
#write.table(reject, file="Results_for_Fig_2.8.csv", sep=",", row.names=F)

# Fig. 2.8
par(mfrow=c(1,1))
plot(glm.Wald ~ b1, data=reject, typ="l", xlab="b1", ylab="Fraction rejected", ylim=c(0,.6), col="green")
lines(glm.LRT ~ b1, data=reject, col="red")
lines(glm.boot0 ~ b1, data=reject, col="black")
lines(lm ~ b1, data=reject, col="blue")
#lines(glm.boot ~ b1, data=reject, col="orange")
lines(logistf ~ b1, data=reject, col="purple")
lines(c(-.1,1.1), c(.05, .05), lty=2)
legend(0,.6,legend=c("Wald", "LRT", "Boot(b1)", "Boot(H0)", "LM", "logistf"), col=c("green","red","orange","black", "blue", "purple"), lty=1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6. As an exercise, produce figure 2.10 showing the power curves for the grouse data assuming that the station-level values of WIND are given by MEAN_WIND for all stations within the same route. This can be done easily by modifying the code provided for making figure 2.9.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
b0 <- -1.0632
b1 <- -0.4064
sd.route <- 1.072

b1.list <- -c(0, .1, .2, .3, .4, .5)

d.sim <- d
w.sim <- w

nsims <- 2000
reject <- data.frame(b1.true=array(0, dim=nsims*length(b1.list)), w.lm=NA, w.glm=NA, w.glm.quasi=NA, w.glmm=NA, d.lm=NA, d.glm=NA, d.glm.anova=NA, d.lmm=NA, d.glmm=NA)

# I set these simulations up a little differently from those in section 2.6. Rather than computing the rejection rate for each b1, here I have kept all of the P-values for all of the simulations and used the aggregate() function to extract P-values from the entire table.

#Get mean_wind values into d.sim 


library(tidyverse)
mean_wind <- d.sim %>%
  group_by(ROUTE) %>%
  summarise(MEAN_WIND = mean(WIND, na.rm = TRUE))

new.d.sim<- left_join(d.sim, mean_wind, by="ROUTE")


i <- 0
for(b1 in b1.list) for(j in 1:nsims){
  i <- i + 1
  reject$b1.true[i] <- b1
  
  d.sim$RUGR <- simulate.d.glmm(d = d, sd = sd.route, b0 = b0, b1 = b1)
  w.sim$RUGR <- aggregate(d.sim$RUGR, by = list(d.sim$ROUTE), FUN = sum)[,2]
  w.sim$SUCCESS <- cbind(w.sim$RUGR, w.sim$STATIONS - w.sim$RUGR)
  
  z.w.lm <- lm(RUGR ~ MEAN_WIND, data=w.sim)
  z.w.glm <- glm(SUCCESS ~ MEAN_WIND, family = binomial, data=w.sim)
  z.w.glm.quasi <- glm(SUCCESS ~ MEAN_WIND, family = quasibinomial, data=w.sim)
  z.w.glmm <- glmer(SUCCESS ~ MEAN_WIND + (1 | ROUTE), family = binomial, data=w.sim, control=glmerControl(calc.derivs=FALSE))
  
  z.d.lm <- lm(RUGR ~ MEAN_WIND, data=new.d.sim)
  z.d.glm <- glm(RUGR ~ MEAN_WIND, family = "binomial", data=new.d.sim)
  z.d.glm.anova <- glm(RUGR ~ MEAN_WIND + ROUTE, family = "binomial", data=new.d.sim)
  z.d.glmm <- lmer(RUGR ~ MEAN_WIND + (1 | ROUTE), data=new.d.sim)
  z.d.lmm <- glmer(RUGR ~ MEAN_WIND + (1 | ROUTE), family = "binomial", data=new.d.sim, control=glmerControl(calc.derivs=FALSE))
  
  reject$w.lm[i] <- summary(z.w.lm)$coef[2,4]	
  reject$w.glm[i] <- summary(z.w.glm)$coef[2,4]	
  reject$w.glm.quasi[i] <- summary(z.w.glm.quasi)$coef[2,4]
  reject$w.glmm[i] <- summary(z.w.glmm)$coef[2,4]
  
  reject$d.lm[i] <- summary(z.d.lm)$coef[2,4]
  reject$d.glm[i] <- summary(z.d.glm)$coef[2,4]
  reject$d.glm.anova[i] <- summary(z.d.glm.anova)$coef[2,4]
  reject$d.glmm[i] <- summary(z.d.glmm)$coef[2,5]
  reject$d.lmm[i] <- summary(z.d.lmm)$coef[2,4]
}
# Because the code takes so long to run, I saved it.
#write.table(reject, file="Results_for_Fig_2.9.csv", sep=",", row.names=F)

# Route-level methods
r <- (reject$w.lm < 0.05)
w.lm <- aggregate(r, by = list(reject$b1.true), FUN = mean)
names(w.lm) <- c("b1", "rejected")

r <- (reject$w.glm < 0.05)
w.glm <- aggregate(r, by = list(reject$b1.true), FUN = mean)
names(w.glm) <- c("b1", "rejected")

r <- (reject$w.glm.quasi < 0.05)
w.glm.quasi <- aggregate(r, by = list(reject$b1.true), FUN = mean)
names(w.glm.quasi) <- c("b1", "rejected")

r <- (reject$w.glmm < 0.05)
w.glmm <- aggregate(r, by = list(reject$b1.true), FUN = mean)
names(w.glmm) <- c("b1", "rejected")

# Fig. 2.9
par(mfrow=c(1,2))
plot(rejected ~ b1, data=w.lm, typ="l", main="Aggregated data", ylab="Fraction rejected", ylim=c(0, 1))
lines(rejected ~ b1, data=w.glm, col="blue")
lines(rejected ~ b1, data=w.glm.quasi, col="green")
lines(rejected ~ b1, data=w.glmm, col="red")
lines(c(-10,0), c(.05,.05), lty=2)
legend(-.25,1,legend=c("LM", "GLM", "GLM.quasi", "GLMM"), col=c("black","blue","green","red"), lty=1)

# Station-level methods
r <- (reject$d.lm < 0.05)
d.lm <- aggregate(r, by = list(reject$b1.true), FUN = mean)
names(d.lm) <- c("b1", "rejected")

r <- (reject$d.glm < 0.05)
d.glm <- aggregate(r, by = list(reject$b1.true), FUN = mean)
names(d.glm) <- c("b1", "rejected")

r <- (reject$d.glm.anova < 0.05)
d.glm.anova <- aggregate(r, by = list(reject$b1.true), FUN = mean)
names(d.glm.anova) <- c("b1", "rejected")

r <- (reject$d.glmm < 0.05)
d.glmm <- aggregate(r, by = list(reject$b1.true), FUN = mean)
names(d.glmm) <- c("b1", "rejected")

r <- (reject$d.lmm < 0.05)
d.lmm <- aggregate(r, by = list(reject$b1.true), FUN = mean)
names(d.lmm) <- c("b1", "rejected")

plot(rejected ~ b1, data=d.lm, typ="l", main="Hierarchical data", ylab="Fraction rejected", ylim=c(0, 1))
lines(rejected ~ b1, data=d.glm, col="blue")
lines(rejected ~ b1, data=d.glm.anova, col="turquoise")
lines(rejected ~ b1, data=d.glmm, col="red")
lines(rejected ~ b1, data=d.lmm, col="orange")
lines(c(-10,0), c(.05,.05), lty=2)
legend(-.25,1,legend=c("LM", "GLM", "GLM.anova", "GLMM", "LMM"), col=c("black","blue","turquoise","red","orange"), lty=1)

#Why are we not rejecting anything? 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 2.7.3 Do hierarchical methods have more power?
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Yes - but not because they have more points, rather they have more information contained within the methods
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################################################################
# Packages
#################################################################

# Packages you will need for Chapter 2:
library(lme4)
library(lmerTest)
# The mvtnorm package has a function for randomly simulating from a multivariate normal distribution.
library(mvtnorm)

# This is the inverse logit function that is used later for simulations
inv.logit <- function(x){
  1/(1 + exp(-x))
}


#################################################################
# 2.3 Estimators
#################################################################

######
# 2.3.1 Maximum Likelihood (ML) estimators 

# This investigates the ML estimator for p, the parameter of a Bernoulli process
p.list <- .01*(1:99)

# Example outcome of the coin flips
x <- c(1,1,0,0,0,0,0,0,0,0)

# Compute the log Likelihood, LL
LL <- data.frame(p=p.list, val=0)
for(i in 1:length(p.list)) {
  LL$p[i] <- p.list[i]
  LL$val[i] <- sum(log(x*p.list[i] + (1-x)*(1-p.list[i])))
}

# Note that this is the same as
m <- 2
LL <- data.frame(p=p.list, val=0)
for(i in 1:length(p.list)) {
  LL$p[i] <- p.list[i]
  LL$val[i] <- m*log(p.list[i]) + (10-m)*log(1-p.list[i])
}

# Fig. 2.1
# Plot the likelihood function LL
plot(val ~ p, data=LL, typ="l", ylab="log L", ylim=c(-35,-5))
lines(c(mean(x), mean(x)), c(-50, max(LL$val)), col=1, lty=2)

#################################################################
# 2.4 Properties of estimators
#################################################################

######
# 2.4.1 Bias 

# This simulates nsims linear regressions and estimates the regression coefficients for each.
n <- 10
b0 <- 1
b1 <- 1
sd.e <- 1

nsims <- 5000
estimate <- data.frame(b0=array(0,dim=nsims), b1=0, b1.P=0, SSRn=0, MSE=0)
for(i in 1:nsims){
  x <- rnorm(n=n, mean=0, sd=1)
  Y <- b0 + b1*x + rnorm(n=n, mean=0, sd=sd.e)
  mod <- lm(Y ~ x)
  estimate$b0[i] <- mod$coef[1]
  estimate$b1[i] <- mod$coef[2]
  estimate$b1.P[i] <- summary(mod)$coef[2,4]
  
  # This gives the total sum-of-squared residuals
  estimate$SSRn[i] <- mean(mod$resid^2)
  
  # This gives the mean squared error (the mean of the sum-of-squared residuals corrected for the degrees of freedom, n/mod$df.residual)
  estimate$MSE[i] <- mean(mod$resid^2)*n/mod$df.residual
}

# Fig. 2.2
par(mfrow = c(1,2))

plot(Y ~ x, main="Single example", ylim=c(b0-3*b1,b0+3*b1))

# The black line gives the fitted line
lines(mod$fitted.values[order(x)] ~ x[order(x)])

# The red line is the line b0 + b1*x with the parameters used to simulation data
lines(b0 + b1*x[order(x)] ~ x[order(x)], col="red")

# This shows the estimator of b1
hist(estimate$b1, xlab="Estimates of b1", main = paste('mean = ', round(mean(estimate$b1),3), sep=''), freq=F, breaks=40)

# This is the formula for plotting the probability density function of a t-distribution that should give the distribution of the estimator of b1
SE <- sd.e*(n/(n-1))*n^-.5
lines(b1 + .01*(-1000:1000), 1/SE*dt(1/SE*.01*(-1000:1000), df=mod$df.residual), col="red")

######
# 2.4.2 Efficiency (Precision) 

# Fig. 2.3
# This plots the ML estimates of the error variance, sd.e^2 (the sum-of-squared residuals divided by n, SSRn)
par(mfrow = c(1,2))
hist(estimate$SSRn, main = paste('mean = ', round(mean(estimate$SSRn),3), sep=''), freq=F, breaks=40, xlim=c(0, max(estimate$SSRn)), xlab="ML estimates")

# This is the formula for plotting the probability density function of a chi-square distribution that approximately gives the distribution of SSRn
lines(sd.e^2*.05*(0:80), (n/sd.e^2)*dchisq(n*.05*(0:80), df=mod$df.residual), col="red")

# This plots the REML estimates of the error variance, sd.e^2 (the Mean Squared Error, MSE)
hist(estimate$MSE, main = paste('mean = ', round(mean(estimate$MSE),3), sep=''), freq=F, breaks=40, xlim=c(0, max(estimate$SSRn)), xlab="REML estimates (MSE)")
lines(sd.e^2*.05*(0:80), ((n-2)/sd.e^2)*dchisq(n*.05*(0:80), df=mod$df.residual), col="red")

# This is the formula for plotting the probability density function of a chi-square distribution that approximately gives the distribution of MSE
lines(sd.e^2*.05*(0:80), ((n-2)/sd.e^2)*dchisq((n-2)*.05*(0:80), df=mod$df.residual), col="green")

legend(x = 1.5, y = .8, c("ML", "REML"), col = c("red","green"), lty = 1)

######
# 2.4.3 Consistency 

# This computes the standard deviation of the estimator of b1 (i.e., the standard error of the estimate of b1) by simulating nsims datasets.
b0 <- 0
b1 <- 1
sd.e <- 1

# this contains the values of the sample size n for different sets of simulations
nrange <- c(10, 20, 50, 100, 200, 500, 1000)

nsims <- 5000
estimate <- data.frame(n=array(0, dim=nsims*length(nrange)), b0=0, b1=0)

counter <- 0
for(n in nrange) for(i in 1:nsims){
  counter <- counter + 1
  x <- rnorm(n=n, mean=0, sd=1)
  Y <- b0 + b1*x + rnorm(n=n, mean=0, sd=sd.e)
  mod <- lm(Y ~ x)
  estimate$n[counter] <- n
  estimate$b0[counter] <- mod$coef[1]
  estimate$b1[counter] <- mod$coef[2]
}
consistency <- aggregate(estimate$b1, by = list(estimate$n), FUN = sd)
names(consistency) <- c("n", "sd.b1")

# Fig. 2.4
par(mfrow=c(1,1))
plot(sd.b1 ~ n, data = consistency, ylab="SD of the estimator",typ="l", ylim=c(.01, max(consistency$sd.b1)))

# This can be used to add lines for additional values of b1
lines(sd.b1 ~ n, data=consistency, col=3)

#################################################################
# 2.5 Hypothesis testing
#################################################################

######
# 2.5.1 Type I errors 

# This counts the proportion of simulations for which the null hypothesis H0:b1=0 is rejected with a value of b1 = 0 in the simulations.

# set b1 = 0 to match the null hypothsis
b0 <- 1
b1 <- 0
sd.e <- 1

# Try this for different sample sizes n
n.list <- c(20,30,50,100)
nsims <- 10000

reject <- data.frame(n=n.list, lm=NA)
n.i <- 0

# loop over the different sample sizes
for(n in n.list){
  n.i <- n.i + 1
  
  # set up a data.frame to collect the P-values for each sample size n
  Pvalues <- data.frame(lm=array(NA,dim=nsims))
  for(i in 1:nsims){
    x <- rnorm(n=n, mean=0, sd=1)
    Y <- b0 + b1*x + rnorm(n=n, mean=0, sd=sd.e)
    mod <- lm(Y ~ x)
    
    Pvalues$lm[i] <- summary(mod)$coef[2,4]
  }
  
  # Calculate the proportion of P-values < 0.05
  reject$lm[n.i] <- mean(Pvalues$lm < 0.05)
}
round(reject, digits=3)

######
# 2.5.2 Power 

# This counts the proportion of simulations for which the null hypothesis H0:b1=0 is rejected when b1 ranges from 0 to 0.5 in the simulations.
b0 <- 1

# start with b1 until the null hypothesis (b1 = 0) and then increase it from the null hypothesis
b1.list <- c(0, .1, .2, .3, .4, .5)
sd.e <- 1

n.list <- c(100, 50, 30, 20)
nsims <- 10000

reject <- data.frame(n=rep(NA, length(n.list)*length(b1.list)), b1=NA, lm=NA)
counter <- 0
for(n in n.list){
  for(b1 in b1.list){
    counter <- counter + 1
    
    # set up a data.frame to collect the P-values for each sample size n
    Pvalues <- data.frame(lm=array(NA,dim=nsims))
    for(i in 1:nsims){
      x <- rnorm(n=n, mean=0, sd=1)
      Y <- b0 + b1*x + rnorm(n=n, mean=0, sd=sd.e)
      mod <- lm(Y ~ x)
      
      Pvalues$lm[i] <- summary(mod)$coef[2,4]
    }
    reject$n[counter] <- n
    reject$b1[counter] <- b1
    reject$lm[counter] <- mean(Pvalues$lm < 0.05)
    show(reject[counter,])
  }
}

# Fig. 2.5
plot(lm ~ b1, data=reject[reject$n == n.list[1],], ylim = c(0,1), typ="l", xlab="b1", ylab="Fraction rejected")
for(i in 2:4) lines(lm ~ b1, data=reject[reject$n == n.list[i],], col=i)
lines(c(-.1,1.1), c(.05, .05), lty=2)
legend(0,1,legend=c("n=100", "    50", "    30", "    20"), col=1:4, lty=1)

#################################################################
# 2.6 P-values for binary data
#################################################################

# The code below is used to produce the table at the top of this section. It uses code that is presented later in this section, so you don't need to run it. And it takes a little time to run. So go ahead and skip down to section 2.6.3

# Get critical values of the LLR for the parametric bootstrap around H0. These are placed in LLR.crit.
b0 <- 1
b1 <- 0
b2 <- 1.5
var.x <- 1
cov.x <- 0
n.list <- c(20,30,50,100)
nboot <- 10000
LLR.crit <- data.frame(n=n.list, glm=NA, glm.converge=NA)
i.n <- 0
for(n in n.list){
  i.n <- i.n + 1
  LLRdeviance <- data.frame(glm=rep(NA, nboot), glm.converge=NA)
  for(i in 1:nboot){
    x <- rmvnorm(n, sigma = matrix(c(var.x, cov.x, cov.x, var.x), nrow=2))
    x1 <- x[,1]
    x2 <- x[,2]
    Y <- rbinom(n=n, size=1, prob=inv.logit(b0 + b1*x1 + b2*x2))
    
    mod.glm <- glm(Y ~ x1 + x2, family = binomial)
    mod.glm0 <- glm(Y ~ x2, family = binomial)
    
    LLRdeviance$glm[i] <- 2*(logLik(mod.glm) - logLik(mod.glm0))[1]
    LLRdeviance$glm.converge[i] <- mod.glm$converge
  }
  LLR.crit$glm[i.n] <- sort(LLRdeviance$glm)[floor(.95*nboot)]
  LLR.converge <- LLRdeviance$glm[LLRdeviance$glm.converge==T]
  LLR.crit$glm.converge[i.n] <- sort(LLR.converge)[floor(.95*length(LLR.converge))]
  show(LLR.crit[i.n, ])
}
round(LLR.crit, digits=2)
# for cov.e=0
# n  glm glm.converge
# 1  20 5.84         4.85
# 2  30 4.61         4.51
# 3  50 4.09         4.09
# 4 100 4.07         4.07

# Compute the rejection rates for the methods except the parametric bootstrap of b1 (which is done in the following code).
b0 <- 1
b1 <- 0
b2 <- 1.5
var.x <- 1
cov.x <- 0
n.list <- c(20,30,50,100)
nsims <- 10000

reject <- data.frame(n=n.list, glm.Wald=NA, glm.Wald.con=NA, glm.LRT=NA, glm.boot0=NA, glm.boot0.con=NA, lm=NA, glm.con=NA)
i.n <- 0
for(n in n.list){
  i.n <- i.n + 1
  Pvalues <- data.frame(glm.Wald=array(NA,dim=nsims), glm.con=NA, glm.LRT=NA, glm.boot0=NA, lm=NA)
  for(i in 1:nsims){
    x <- rmvnorm(n, sigma = matrix(c(var.x, cov.x, cov.x, var.x), nrow=2))
    x1 <- x[,1]
    x2 <- x[,2]
    Y <- rbinom(n=n, size=1, prob=inv.logit(b0 + b1*x1 + b2*x2))
    
    mod.glm <- glm(Y ~ x1 + x2, family = binomial)
    mod.glm0 <- glm(Y ~ x2, family = binomial)
    mod.lm <- lm(Y ~ x1 + x2)
    
    Pvalues$glm.Wald[i] <- summary(mod.glm)$coef[2,4]
    Pvalues$glm.LRT[i] <- pchisq(2*(logLik(mod.glm) - logLik(mod.glm0)), df=1, lower.tail=F)
    Pvalues$glm.deviance[i] <- 2*(logLik(mod.glm) - logLik(mod.glm0))
    Pvalues$lm[i] <- summary(mod.lm)$coef[2,4]
    
    Pvalues$glm.con[i] <- mod.glm$converge
  }
  reject$glm.Wald[i.n] <- mean(Pvalues$glm.Wald < 0.05)
  reject$lm[i.n] <- mean(Pvalues$lm < 0.05)
  reject$glm.Wald.con[i.n] <- mean(Pvalues$glm.Wald[Pvalues$glm.con==T] < 0.05)
  reject$glm.LRT[i.n] <- mean(Pvalues$glm.LRT < 0.05)
  reject$glm.boot0[i.n] <- mean(Pvalues$glm.deviance > LLR.crit[n.i,2])
  reject$glm.boot0.con[i.n] <- mean(Pvalues$glm.deviance[Pvalues$glm.con==T] > LLR.crit[LLR.crit$n == n,3])
  reject$glm.con[i.n] <- mean(Pvalues$glm.con)
  show(reject[i.n,])
}
round(reject, digits=3)
# n glm.Wald glm.Wald.con glm.LRT glm.boot0 glm.boot0.con    lm glm.con
# 1  20    0.004        0.005   0.103     0.047         0.048 0.050   0.952
# 2  30    0.025        0.026   0.078     0.054         0.054 0.053   0.994
# 3  50    0.042        0.042   0.063     0.055         0.055 0.054   1.000
# 4 100    0.046        0.046   0.055     0.049         0.049 0.049   1.000

# Parametric bootstrap of b1. Note that this takes a lot of time, so don't run it.
b0 <- 1
b1 <- 0
b2 <- 1.5
var.x <- 1
cov.x <- 0
n.list <- c(20,30,50,100)
nsims <- 500
nboot <- 2000

reject <- data.frame(n=n.list, glm.boot=NA)
i.n <- 0
for(n in n.list){
  i.n <- i.n + 1
  Pvalues <- data.frame(glm.boot=rep(NA, nsims))
  for(i in 1:nsims){
    x <- rmvnorm(n, sigma = matrix(c(var.x, cov.x, cov.x, var.x), nrow=2))
    x1 <- x[,1]
    x2 <- x[,2]
    Y <- rbinom(n=n, size=1, prob=inv.logit(b0 + b1*x1 + b2*x2))
    dat <- data.frame(x1=x1, x2=x2, Y=Y)
    
    mod.real <- glm(Y ~ x1 + x2, family = binomial, dat)
    # Parametric bootstrap
    estimate <- data.frame(b1=array(NA,dim=nboot))
    dat.boot <- dat
    for(j in 1:nboot){
      dat.boot$Y <- simulate(mod.real)[[1]]
      mod.boot <- glm(Y ~ x1 + x2, family = binomial, data=dat.boot)
      if(mod.boot$converge == T) estimate$b1[j] <- mod.boot$coef[2]
    }
    Pvalues$glm.boot[i] <- 2*min(mean(estimate$b1 < 0, na.rm=T), mean(estimate$b1 > 0, na.rm=T))
  }
  reject$glm.boot[i.n] <- mean(Pvalues$glm.boot < 0.05, na.rm=T)
  show(reject[i.n,])
}
round(reject, digits=3)
# n  glm.boot
# 1  20 0.1016598
# 2  30 0.0920000
# 3  50 0.0620000
# 4 100 0.0740000

######
#2.6.3 Parametric bootstrap of b1
n <- 30
b0 <- 1
b1 <- 1
b2 <- 0
var.x <- 1
cov.x <- 0
nboot <- 1000

# Simulate the "real" dataset dat. The rmvnorm() function is a multivariate normal, so that but x1 and x2 are normally distributed with variances var.x and covariance cov.x between them.
x <- rmvnorm(n, sigma = matrix(c(var.x, cov.x, cov.x, var.x), nrow=2))
x1 <- x[,1]
x2 <- x[,2]
Y <- rbinom(n=n, size=1, prob=inv.logit(b0 + b1*x1 + b2*x2))
dat <- data.frame(x1=x1, x2=x2, Y=Y)

# Fit dat with the full model.
mod.dat <- glm(Y ~ x1 + x2, family = binomial, data=dat)
summary(mod.dat)

# Perform the parametric bootstrap by simulating data with the full model.
boot <- data.frame(b1=array(NA, dim=nboot), converge=NA)
dat.boot <- dat
for(i in 1:nboot){
  dat.boot$Y <- simulate(mod.dat)[[1]]
  mod.boot <- glm(Y ~ x1 + x2, family = binomial, data=dat.boot)
  boot$converge[i] <- mod.boot$converge
  boot$b1[i] <- mod.boot$coef[2]
}
# Remove the cases when glm() did not converge
boot <- boot[boot$converge == T,]
pvalue <- 2*min(mean(boot$b1 < 0), mean(boot$b1 > 0))
pvalue

######
#2.6.4 Parametric bootstrap of H0

# Fit full and reduced models to the data.
mod.dat <- glm(Y ~ x1 + x2, family = binomial, data=dat)
mod.reduced <- glm(Y ~ x2, family = binomial, data=dat)
dev.dat <- 2*(logLik(mod.dat) - logLik(mod.reduced))[1]

# Simulate the data using the reduced model.
boot0 <- data.frame(LLR=rep(NA, nboot), converge=NA)
dat.boot <- dat
for(i in 1:nboot){
  dat.boot$Y <- simulate(mod.reduced)[[1]]
  mod.glm <- update(mod.dat, data=dat.boot)
  mod.glm0 <- update(mod.reduced, data=dat.boot)
  boot0$dev[i] <- 2*(logLik(mod.glm) - logLik(mod.glm0))[1]
  boot0$converge[i] <- mod.glm$converge
}
# Remove the cases when glm() did not converge
boot0 <- boot0[boot0$converge == T,]
pvalue <- mean(boot0$dev > dev.dat)
pvalue

# This is the p-value from a LRT
pchisq(dev.dat, df=1, lower.tail=F)

# Fig. 2.6 
# Comparing both bootstraps
par(mfrow=c(1,2))

# For small n, it might be necessary to remove poorly fit simulations
boot <- boot[boot$b1 < (mean(boot$b1) + 2.5*sd(boot$b1)),]

# Plot of parametric bootstrap
hist(boot$b1, xlab="Bootstrap of b1", main = paste("b1 = ", round(mod.dat$coef[2],2), ",  Estimate mean = ", round(mean(boot$b1),2), sep=''), freq=F, breaks=40)
lines(c(mod.dat$coef[2],mod.dat$coef[2]), c(0,2), col="red")
lines(c(0,0), c(0,2), col="green")
lines(.1*(-100:100), dnorm(.1*(-100:100), mean=mod.dat$coef[2], sd=summary(mod.dat)$coef[2,2]), col="blue")

# Plot of parametric bootstrap around H0
hist(boot0$dev, xlab="Bootstrap of deviance", main = paste("Deviance = ", round(dev.dat,2), sep=''), freq=F, breaks=40)
lines(c(dev.dat,dev.dat), c(0,1), col="red")
lines(.1*(1:100), dchisq(.1*(1:100), df=1), col="blue")

######
#2.6.5 Why did the LM do so well?

# Simulate the "real" dataset dat with many points under H0:b1=0 but with b2 > 0 and covariance between x1 and x2
n <- 1000
b0 <- 1
b1 <- 0
b2 <- 1.5
var.x <- 1
cov.x <- .8

x <- rmvnorm(n, sigma = matrix(c(var.x, cov.x, cov.x, var.x), nrow=2))
x1 <- x[,1]
x2 <- x[,2]
Y <- rbinom(n=n, size=1, prob=inv.logit(b0 + b1*x1 + b2*x2))
dat <- data.frame(x1=x1, x2=x2, Y=Y)

# Fit the reduced LM model
mod.reduced <- lm(Y ~ x2, data=dat)
dat$r <- mod.reduced$resid

# Fig. 2.7
# Plot the residuals
par(mfrow=c(1,2))
plot(r ~ x1, data=dat, xlab="Values of x1", ylab="Residuals")

# Bin and plot the variances of the residuals
dat <- dat[order(dat$x1),]
bins <- round(n*(0:10)/10)
resid <- data.frame(x1=bins[-1], var=NA)
for(i in 2:11){
  resid$x1[i-1] <- mean(dat$x1[bins[i-1]:bins[i]])
  resid$var[i-1] <- var(dat$r[bins[i-1]:bins[i]])
}
plot(var ~ x1, data=resid, ylim=c(0,max(var)), xlab="Binned values of x1", ylab="Residual variance")

######
#2.6.6 Power in fitting the binary data

# Compute the rejection rates for the methods except the parametric bootstrap of b1 (which is done in the following code). Note that you will need to calculate the distribution of LLR for bootstrap of H0 using the code at the top of section "2.6 P-values for binary data" above.
n <- 30
b0 <- 1
b1 <- 0
b2 <- 1.5
var.x <- 1
cov.x <- 0
b1.list <- c(0, .2, .4, .6, .8, 1)
nsims <- 10000

reject <- data.frame(b1=b1.list, glm.Wald=NA, glm.Wald.con=NA, glm.LRT=NA, glm.boot0=NA, glm.boot0.con=NA, lm=NA, glm.con=NA)
i.b1 <- 0
for(b1 in b1.list){
  i.b1 <- i.b1 + 1
  Pvalues <- data.frame(glm.Wald=array(NA,dim=nsims), glm.con=NA, glm.LRT=NA, glm.boot0=NA, lm=NA)
  for(i in 1:nsims){
    x <- rmvnorm(n, sigma = matrix(c(var.x, cov.x, cov.x, var.x), nrow=2))
    x1 <- x[,1]
    x2 <- x[,2]
    Y <- rbinom(n=n, size=1, prob=inv.logit(b0 + b1*x1 + b2*x2))
    
    mod.glm <- glm(Y ~ x1 + x2, family = binomial)
    mod.glm0 <- glm(Y ~ x2, family = binomial)
    mod.lm <- lm(Y ~ x1 + x2)
    
    Pvalues$glm.Wald[i] <- summary(mod.glm)$coef[2,4]
    Pvalues$glm.LRT[i] <- pchisq(2*(logLik(mod.glm) - logLik(mod.glm0)), df=1, lower.tail=F)
    Pvalues$glm.deviance[i] <- 2*(logLik(mod.glm) - logLik(mod.glm0))
    Pvalues$lm[i] <- summary(mod.lm)$coef[2,4]
    
    Pvalues$glm.con[i] <- mod.glm$converge
  }
  reject$glm.Wald[i.b1] <- mean(Pvalues$glm.Wald < 0.05)
  reject$lm[i.b1] <- mean(Pvalues$lm < 0.05)
  reject$glm.Wald.con[i.b1] <- mean(Pvalues$glm.Wald[Pvalues$glm.con==T] < 0.05)
  reject$glm.LRT[i.b1] <- mean(Pvalues$glm.LRT < 0.05)
  reject$glm.boot0[i.b1] <- mean(Pvalues$glm.deviance > LLR.crit[LLR.crit$n == n,3])
  reject$glm.boot0.con[i.b1] <- mean(Pvalues$glm.deviance[Pvalues$glm.con==T] > LLR.crit[LLR.crit$n == n,3])
  reject$glm.con[i.b1] <- mean(Pvalues$glm.con)
  show(reject[i.b1,])
}
round(reject, digits=3)
# b1 glm.Wald glm.Wald.con glm.LRT glm.boot0 glm.boot0.con     lm glm.con
# 1 0.0   0.0228   0.02290997  0.0755    0.0536    0.05144695 0.0505  0.9952
# 2 0.2   0.0389   0.03912694  0.1001    0.0734    0.07071012 0.0723  0.9942
# 3 0.4   0.0715   0.07190989  0.1620    0.1238    0.12038620 0.1244  0.9943
# 4 0.6   0.1406   0.14160540  0.2748    0.2236    0.21945815 0.2234  0.9929
# 5 0.8   0.2431   0.24478904  0.4064    0.3407    0.33732756 0.3427  0.9931
# 6 1.0   0.3526   0.35576632  0.5500    0.4819    0.47845828 0.4820  0.9911

# Parametric bootstrap of b1. Note that this takes a lot of time to run.
nsims <- 500
nboot <- 2000

reject.boot <- data.frame(b1=b1.list, glm.boot=NA)
i.b1 <- 0
for(b1 in b1.list){
  i.b1 <- i.b1 + 1
  Pvalues <- data.frame(glm.boot=rep(NA, nsims))
  for(i in 1:nsims){
    x <- rmvnorm(n, sigma = matrix(c(var.x, cov.x, cov.x, var.x), nrow=2))
    x1 <- x[,1]
    x2 <- x[,2]
    Y <- rbinom(n=n, size=1, prob=inv.logit(b0 + b1*x1 + b2*x2))
    dat <- data.frame(x1=x1, x2=x2, Y=Y)
    
    mod.real <- glm(Y ~ x1 + x2, family = binomial, dat)
    # Parametric bootstrap
    estimate <- data.frame(b1=array(NA,dim=nboot))
    dat.boot <- dat
    for(j in 1:nboot){
      dat.boot$Y <- simulate(mod.real)[[1]]
      mod.boot <- glm(Y ~ x1 + x2, family = binomial, data=dat.boot)
      if(mod.boot$converge == T) estimate$b1[j] <- mod.boot$coef[2]
    }
    Pvalues$glm.boot[i] <- 2*min(mean(estimate$b1 < 0, na.rm=T), mean(estimate$b1 > 0, na.rm=T))
  }
  reject.boot$glm.boot[i.b1] <- mean(Pvalues$glm.boot < 0.05, na.rm=T)
  show(reject.boot[i.b1,])
}
# b1  glm.boot
# 1 0.0 0.1022044
# 2 0.2 0.1064257
# 3 0.4 0.1943888
# 4 0.6 0.2897384
# 5 0.8 0.4567404
# 6 1.0 0.5688259

reject$glm.boot <- reject.boot$glm.boot
# Because the code takes so long to run, I saved it.
write.table(reject, file="Results_for_Fig_2.8.csv", sep=",", row.names=F)

# Fig. 2.8
par(mfrow=c(1,1))
plot(glm.Wald ~ b1, data=reject, typ="l", xlab="b1", ylab="Fraction rejected", ylim=c(0,.6), col="green")
lines(glm.LRT ~ b1, data=reject, col="red")
lines(glm.boot0 ~ b1, data=reject, col="black")
lines(lm ~ b1, data=reject, col="blue")
lines(glm.boot ~ b1, data=reject, col="orange")
lines(c(-.1,1.1), c(.05, .05), lty=2)
legend(0,.6,legend=c("Wald", "LRT", "Boot(b1)", "Boot(H0)", "LM"), col=c("green","red","orange","black", "blue"), lty=1)


#################################################################
# 2.7 Example data: grouse
#################################################################

# read data
d <- read.csv(file="grouse_data.csv", header=T)

# STATION and ROUTE were uploaded as integers; this converts them to factors.
d$STATION <- as.factor(d$STATION)
d$ROUTE <- as.factor(d$ROUTE)

# Combine all values for each ROUTE.
w <- data.frame(aggregate(cbind(d$GROUSE, d$WIND, d$LAT, d$LONG), by = list(d$ROUTE), FUN = mean))

# For clarity, I've added the column names.
names(w) <- c("ROUTE","MEAN_GROUSE","MEAN_WIND","LAT","LONG")

# Finally, I want the count of the number of STATIONS per ROUTE, and the number of GROUSE.
w$GROUSE <- aggregate(d$GROUSE, by = list(d$ROUTE), FUN = sum)[,2]
w$STATIONS <- aggregate(array(1, c(nrow(d), 1)), by = list(d$ROUTE), FUN = sum)[,2]

######
# 2.7.1 Simulating the grouse data

# This function simulates grouse data given an input dataset d (which is used to give the structure of the dataset and values of WIND), and values of parameters b0, b1, and sd (residual standard deviation).
simulate.d.glmm <- function (d, b0, b1, sd) {
  # This is the inverse logit function
  inv.logit <- function(x){
    1/(1 + exp(-x))
  }
  d.sim.Y <- array(-1, dim=dim(d)[1])
  for(i in levels(d$ROUTE)){
    dd <- d[d$ROUTE == i,]
    nn <- dim(dd)[1]
    Z <- b0 + b1*dd$WIND + rnorm(n=1, mean=0, sd=sd)
    p <- inv.logit(Z)
    d.sim.Y[d$ROUTE == i] <- rbinom(n=nn, size=1, prob=p)
  }
  return(d.sim.Y)
}

######
# 2.7.2 Power curves for the grouse data
# Warning: this can take a long time (>1 hour) to run. To just check out the code, you can decrease nsims to 20.

# These are parameter values estimated from the data.
b0 <- -1.0632
b1 <- -0.4064
sd.route <- 1.072

b1.list <- -c(0, .1, .2, .3, .4, .5)

d.sim <- d
w.sim <- w

nsims <- 2000
reject <- data.frame(b1.true=array(0, dim=nsims*length(b1.list)), w.lm=NA, w.glm=NA, w.glm.quasi=NA, w.glmm=NA, d.lm=NA, d.glm=NA, d.glm.anova=NA, d.lmm=NA, d.glmm=NA)

# I set these simulations up a little differently from those in section 2.6. Rather than computing the rejection rate for each b1, here I have kept all of the P-values for all of the simulations and used the aggregate() function to extract P-values from the entire table.
i <- 0
for(b1 in b1.list) for(j in 1:nsims){
  i <- i + 1
  reject$b1.true[i] <- b1
  
  d.sim$RUGR <- simulate.d.glmm(d = d, sd = sd.route, b0 = b0, b1 = b1)
  w.sim$RUGR <- aggregate(d.sim$RUGR, by = list(d.sim$ROUTE), FUN = sum)[,2]
  w.sim$SUCCESS <- cbind(w.sim$RUGR, w.sim$STATIONS - w.sim$RUGR)
  
  z.w.lm <- lm(RUGR ~ MEAN_WIND, data=w.sim)
  z.w.glm <- glm(SUCCESS ~ MEAN_WIND, family = binomial, data=w.sim)
  z.w.glm.quasi <- glm(SUCCESS ~ MEAN_WIND, family = quasibinomial, data=w.sim)
  z.w.glmm <- glmer(SUCCESS ~ MEAN_WIND + (1 | ROUTE), family = binomial, data=w.sim, control=glmerControl(calc.derivs=FALSE))
  
  z.d.lm <- lm(RUGR ~ WIND, data=d.sim)
  z.d.glm <- glm(RUGR ~ WIND, family = "binomial", data=d.sim)
  z.d.glm.anova <- glm(RUGR ~ WIND + ROUTE, family = "binomial", data=d.sim)
  z.d.glmm <- lmer(RUGR ~ WIND + (1 | ROUTE), data=d.sim)
  z.d.lmm <- glmer(RUGR ~ WIND + (1 | ROUTE), family = "binomial", data=d.sim, control=glmerControl(calc.derivs=FALSE))
  
  reject$w.lm[i] <- summary(z.w.lm)$coef[2,4]	
  reject$w.glm[i] <- summary(z.w.glm)$coef[2,4]	
  reject$w.glm.quasi[i] <- summary(z.w.glm.quasi)$coef[2,4]
  reject$w.glmm[i] <- summary(z.w.glmm)$coef[2,4]
  
  reject$d.lm[i] <- summary(z.d.lm)$coef[2,4]
  reject$d.glm[i] <- summary(z.d.glm)$coef[2,4]
  reject$d.glm.anova[i] <- summary(z.d.glm.anova)$coef[2,4]
  reject$d.glmm[i] <- summary(z.d.glmm)$coef[2,5]
  reject$d.lmm[i] <- summary(z.d.lmm)$coef[2,4]
}
# Because the code takes so long to run, I saved it.
write.table(reject, file="Results_for_Fig_2.9.csv", sep=",", row.names=F)

# Route-level methods
r <- (reject$w.lm < 0.05)
w.lm <- aggregate(r, by = list(reject$b1.true), FUN = mean)
names(w.lm) <- c("b1", "rejected")

r <- (reject$w.glm < 0.05)
w.glm <- aggregate(r, by = list(reject$b1.true), FUN = mean)
names(w.glm) <- c("b1", "rejected")

r <- (reject$w.glm.quasi < 0.05)
w.glm.quasi <- aggregate(r, by = list(reject$b1.true), FUN = mean)
names(w.glm.quasi) <- c("b1", "rejected")

r <- (reject$w.glmm < 0.05)
w.glmm <- aggregate(r, by = list(reject$b1.true), FUN = mean)
names(w.glmm) <- c("b1", "rejected")

# Fig. 2.9
par(mfrow=c(1,2))
plot(rejected ~ b1, data=w.lm, typ="l", main="Aggregated data", ylab="Fraction rejected", ylim=c(0, 1))
lines(rejected ~ b1, data=w.glm, col="blue")
lines(rejected ~ b1, data=w.glm.quasi, col="green")
lines(rejected ~ b1, data=w.glmm, col="red")
lines(c(-10,0), c(.05,.05), lty=2)
legend(-.25,1,legend=c("LM", "GLM", "GLM.quasi", "GLMM"), col=c("black","blue","green","red"), lty=1)

# Station-level methods
r <- (reject$d.lm < 0.05)
d.lm <- aggregate(r, by = list(reject$b1.true), FUN = mean)
names(d.lm) <- c("b1", "rejected")

r <- (reject$d.glm < 0.05)
d.glm <- aggregate(r, by = list(reject$b1.true), FUN = mean)
names(d.glm) <- c("b1", "rejected")

r <- (reject$d.glm.anova < 0.05)
d.glm.anova <- aggregate(r, by = list(reject$b1.true), FUN = mean)
names(d.glm.anova) <- c("b1", "rejected")

r <- (reject$d.glmm < 0.05)
d.glmm <- aggregate(r, by = list(reject$b1.true), FUN = mean)
names(d.glmm) <- c("b1", "rejected")

r <- (reject$d.lmm < 0.05)
d.lmm <- aggregate(r, by = list(reject$b1.true), FUN = mean)
names(d.lmm) <- c("b1", "rejected")

plot(rejected ~ b1, data=d.lm, typ="l", main="Hierarchical data", ylab="Fraction rejected", ylim=c(0, 1))
lines(rejected ~ b1, data=d.glm, col="blue")
lines(rejected ~ b1, data=d.glm.anova, col="turquoise")
lines(rejected ~ b1, data=d.glmm, col="red")
lines(rejected ~ b1, data=d.lmm, col="orange")
lines(c(-10,0), c(.05,.05), lty=2)
legend(-.25,1,legend=c("LM", "GLM", "GLM.anova", "GLMM", "LMM"), col=c("black","blue","turquoise","red","orange"), lty=1)


