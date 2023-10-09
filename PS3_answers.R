# PS3
# due 3 October, 2023 by midnight

# This problem set is to check to make sure you are up and running in R, and it begins to explore the ruffed grouse data. It is based on correlated_data Ch1 ("Ch_1_10Aug18.R") sections 1.5 to 1.7.

# For the homework, EACH INDIVIDUAL SHOULD TURN IN A SEPARATE .R FILE. Start with the code from correlated_data ("Ch_1_10Aug18.R") and add to it anything you need. Identify you new code (so I can find it) by placing it between marker rows #~~~~~~~~~~~~~~~~~~~~~~~~~~~. I will get your full answers to the questions by randomly asking people in class, so there is no need for you to write them down.

# 1. What questions do you have about the material in section 1.5 to 1.7? What needs more explanation? I'm serious about asking this question, because I want to improve the book. (NOTE: This requires NO NEW R CODE.)

# 2. Suppose Y1, Y2, and Y3 are independent normal random variables with variances 1, 2, and 3, respectively. You can think of (Y1 + Y3) and (Y2 + Y3) as the number of grouse observed in station 1 and station 2 of route 1: Y3 gives the random effect that is common to all stations within route 1, while Y1 and Y2 give the station-specific random effects for stations 1 and 2 within route 1. What is the covariance between (Y1 + Y3) and (Y2 + Y3)?  If you want, check your answer with a simulation.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
predicted.cov <- 3


Y1=rnorm(10000, mean=0, sd=1)
Y2=rnorm(10000, mean=0, sd=sqrt(2))
Y3=rnorm(10000, mean=0, sd=sqrt(3))
Y4=Y1+Y3
Y5=Y2+Y3
cov(Y4, Y5)
predicted.cov

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Further exploration of covariance matrices

# Suppose you have 3 groups (routes) with 4 samples (stations) per group.
ngroups <- 3
nsamples.per.group <- 4
nsamples <- ngroups * nsamples.per.group

# You want to estimate the covariance matrix for the 12 samples. Suppose the outcome (number of grouse) is given by

#	Z <- b0 + beta.route + e

# where b0 is the mean number of grouse per station, beta.route is route-specific difference from the overall mean (b0) for the 3 groups, and e is a randome variable giving the sample(station)-specific variation. 

# beta.route is assumed to be normally distributed with mean zero and standard deviation sd.beta.route

# set b0 to zero and allow variation among routes using var.beta > 0.
b0 <- 0
var.beta.route <- 3
var.e <- 1

# create a group variable
group <- as.factor(rep(1:ngroups, each=nsamples.per.group))

# perform the simulations nrep times
nrep <- 100000

# create a matrix to collect the simulation results
sim.results <- matrix(NA, nrow = nrep, ncol = nsamples)

# loops of nreps to collect values of Z
for(i.rep in 1:nrep){
  e <- rnorm(nsamples, mean = 0, sd = var.e^.5)
  beta <- rnorm(ngroups, mean = 0, sd = var.beta.route^.5)
  beta.route <- rep(beta, each = nsamples.per.group)
  
  # simulation
  Z <- b0 + beta.route + e
  
  # place Z into the rows of sim.results
  sim.results[i.rep,] <- Z
}
# This is what the simulations look like
head(sim.results)

# This computes the covariance of all the values from sim.results, i.e., the covariance matrix. The round() function makes things easier to see.
round(cov(sim.results), digits = 1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# code from Ch_1_10Aug18.R

#################################################################
# 1.3 Dataset
#################################################################

# Metadata for "grouse_data.csv"

# ROUTE
# IDs for 50 roadside routes. These data are simulated to have similar characteristics as the original, real data.

# STATION
# IDs for each survey station, with up to 8 STATIONS per ROUTE.

# LAT
# X coordinate of survey station. UTM zone 15N NAD83 datum.

# LONG
# Y coordinate of survey station. UTM zone 15N NAD83 datum.

# WIND
# Wind speed (km/hour) recorded at 1.4m above ground at the end of the survey.

# TEMP
# Temperature (°C) recorded at 1.4m above ground at the end of the survey.

# GROUSE
# Detection/non-detection of Ruffed Grouse (1 = detected, 0 = not detected).

# read data
d <- read.csv(file="grouse_data.csv", header=T)

# STATION and ROUTE were uploaded as integers; this converts them to factors.
d$STATION <- as.factor(d$STATION)
d$ROUTE <- as.factor(d$ROUTE)

# To see the first 20 rows of data, you can use:
head(d, 20)

# Combine all values for each ROUTE.
w <- data.frame(aggregate(cbind(d$GROUSE, d$WIND, d$LAT, d$LONG), by = list(d$ROUTE), FUN = mean))

# For clarity, I've added the column names.
names(w) <- c("ROUTE","MEAN_GROUSE","MEAN_WIND","LAT","LONG")

# Finally, I want the count of the number of STATIONS per ROUTE, and the number of GROUSE.
w$GROUSE <- aggregate(d$GROUSE, by = list(d$ROUTE), FUN = sum)[,2]
w$STATIONS <- aggregate(array(1, c(nrow(d), 1)), by = list(d$ROUTE), FUN = sum)[,2]

# This is what the aggregated data look like:
head(w)

# Fig. 1.1
# To take a look at the data, this plots MEAN_GROUSE against LAT and LONG, using the size of the point to give the mean number of observations per route. "cex" gives the point size, and I've added 0.1 to w$MEAN_GROUSE so that routes with zero counts are still shown
par(mfrow=c(1,3))
plot(LAT ~ LONG, data=w, pch=1, cex=5*(MEAN_GROUSE+.1), xaxt="n", yaxt="n", xlab="Longitude", ylab="Latitude")

# This histogram shows the variability in grouse counts per ROUTE.
hist(w$MEAN_GROUSE, main="", xlab="MEAN_GROUSE per ROUTE")

# It looks like there are fewer GROUSE seen when the wind speed is higher
plot(MEAN_GROUSE ~ MEAN_WIND, data=w)

#########################################################
# 1.5 Station-level analyses
#########################################################
library(lme4)

######
# 1.5.1 LM at the station level
summary(lm(GROUSE ~ WIND, data=d))

######
# 1.5.2 GLM at the station level
summary(glm(GROUSE ~ WIND, data=d, family=binomial))

# quasi-GLM: This is not appropriate, because for binary (0 or 1) data, the variance has to be m(1-m)

######
# 1.5.3 GLM with ROUTE as a factor
library(car)
Anova(glm(GROUSE ~ WIND + ROUTE, data=d, family=binomial))

######
# 1.5.4 GLMM: 
summary(glmer(GROUSE ~ WIND + (1 | ROUTE), data=d, family=binomial))

# Demonstration of the structure of a mixed model by comparing the glm having ROUTEs as factors with the glmm having ROUTEs as a random variable.
z.glm <- glm(GROUSE ~ 0 + WIND + ROUTE, data=d, family=binomial)

z.glmm <- glmer(GROUSE ~ WIND + (1 | ROUTE), data=d, family=binomial)

# Fig. 1.3
par(mfrow=c(1,2), mai=c(1,.8,.5,.1), mgp=c(3,.5,0))

# This plots the 50 values of route-factor coefficients
hist(z.glm$coef[2:51], breaks=(-20:20), freq=F,  main="Estimates from GLM", xlab="b1")

# This plots the 50 values of the random effects associated with the 50 routes
hist(coef(z.glmm)$ROUTE[[1]], breaks=(-20:20), freq=F, main="Estimates from GLMM", xlab="b1")
lines(.02*(-1000:1000), dnorm(.02*(-1000:1000), mean=fixef(z.glmm)['WIND'], sd=summary(z.glmm)$var[[1]][1]^.5), col="red")

# This is the variance of the ROUTE coefficients that can be compared to the ROUTE random effect in the GLMM. It is much larger.
var(z.glm$coef[2:51])
# This is due to the outlier estimates of the ROUTE coefficents of zero

# Produce a data.frame d.no0 that only contains ROUTES with at least one grouse observation, and then enalyze it with glm() and glmer(). The results for the glm() and glmm() are much more similar.
w.no0 <- w[w$MEAN_GROUSE > 0,]
d.no0 <- d[is.element(d$ROUTE, w.no0$ROUTE),]

z.glm.no0 <- glm(GROUSE ~ WIND + ROUTE, data=d.no0, family=binomial)
Anova(z.glm.no0)
var(z.glm.no0$coef[2:length(z.glm.no0$coef)])

# Note that the variance in fixed effects for routes in the GLM is now much closer to (but still larger than) the variance of the random effets term for ROUTE in the GLMM.
summary(glmer(GROUSE ~ WIND + (1 | ROUTE), data=d.no0, family=binomial))

######
# 1.5.5 LMM: Even though this is a linear model, it incorporates the non-independence of STATIONS within ROUTES.
summary(lmer(GROUSE ~ WIND + (1 | ROUTE), data=d))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3.  (from book) Fit a LM to the station-level data treating route as a factor. This is similar to the LMM treating route as a random effect (subsection 1.5.5). Do they give similar results about the effect of WIND on GROUSE? Is the comparison between these LM and LMM similar to the comparison between the GLM treating route as a factor (subsection 1.5.3) and the GLMM treating route as a random effect (subsection 1.5.4)? Produce a figure like figure 1.3. Does this explain your answer in the same way figure 1.3 explains the difference between the GLM and GLMM?

summary(lmer(GROUSE ~ WIND + (1 | ROUTE), data=d))
library(lmerTest)
summary(lmer(GROUSE ~ WIND + (1 | ROUTE), data=d))

summary(lm(GROUSE ~ WIND + ROUTE, data=d))
Anova(lm(GROUSE ~ WIND + ROUTE, data=d))

# Note that the P-value for WIND is higher for the LM, indicating a loss of power relative to the LMM.

# To get the equivalent of random effects for ROUTE using the LM, you need to remove the intercept using "0 + ". This gives you the value for all routes, the the "residual" route effect is these values minus their mean value (used in the hist() plot).
z.lm <- lm(GROUSE ~ 0 + WIND + ROUTE, data=d)
z.lmm <- lmer(GROUSE ~ 1 + WIND + (1 | ROUTE), data=d)

# Fig. similar to 1.3. 
par(mfrow=c(1,2), mai=c(1,.8,.5,.1), mgp=c(3,.5,0))
hist(z.lm$coef[2:51] - mean(z.lm$coef[2:51]), breaks=.1*(-20:20), freq=F,  main="Estimates from LM", xlab="b1 (residual route effects)")
hist(unlist(ranef(z.lmm)), breaks=.1*(-20:20), freq=F, main="Estimates from LMM", xlab="b1 (random effects for route)")
lines(.02*(-1000:1000), dnorm(.02*(-1000:1000), mean=fixef(z.lmm)['WIND'], sd=summary(z.lmm)$var[[1]][1]^.5), col="red")

# There is more variation in the predicted route coefficients for the LM than the LMM. 
var(z.lm$coef[2:51])
var(unlist(ranef(z.lmm)))

# The patterns in this figure aren't like figure 1.3; in both the LM and LMM, the values of b1 are clustered around zero. The loss of power to detect the effect of WIND in the LM is likely due to the number of parameters that must be estimated. Even though this leads to a better model fit, the calculation of the P-value penalizes all of the parameters that are used in the model.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. In section 1.4.5 there is a simulation for the GLMM model of route-level data. How could you use a model (it has to be a different model) for station-level data to test whether the P-values given by the LMM (section 1.5.5) are correct? This question is similar to Question 4 in Problem Set 2.

# This is the inverse logit function
inv.logit <- function(x){
  1/(1 + exp(-x))
}

# Set up parameters for the simulation. There are 50 groups with 8 samples per group.
ngroups <- 50
nsamples.per.group <- 8
nsamples <- ngroups * nsamples.per.group

# set the coefficients to zero and allow variation among routes using sd.beta > 0.
b0 <- 0
b1 <- 0
sd.beta <- 1

# create a group variable
group <- as.factor(rep(1:ngroups, each=nsamples.per.group))

# perform the simulations nrep times
nrep <- 2000
test<- data.frame(rep=1:nrep)
for(rep in 1:nrep){
  x <- rnorm(nsamples, mean=0, sd=1)
  e <- rnorm(nsamples, mean=0, sd=1)
  beta <- rnorm(ngroups, mean=0, sd=sd.beta)
  x.beta <- rep(beta, each=nsamples.per.group)
  
  # simulation
  Z <- b0 + x * b1 + x.beta + e
  
  p <- inv.logit(Z)
  Y <- rbinom(nsamples, size=1, prob=p)
  
  # These use the summary() function for extracting P-values
  test$p.value.lmer[rep] <- summary(lmer(Y ~ x + (1|group)))$coef[2,5]
  test$p.value.glmer[rep] <- summary(glmer(Y ~ x + (1|group), family=binomial))$coef[2,4]
}

# Histograms of the P-values. Values for the LMM are in gray and for the GLMM in red.
par(mfrow=c(1,1))
hist(test$p.value.lmer, breaks=20, main=paste0("lmer: ", mean(test$p.value.lmer < 0.05), " glmer: ", mean(test$p.value.glmer < 0.05)))
hist(test$p.value.glmer, breaks=20, add=T, col="red")


#~~~~~~~~~~~~~~~~~
# The code can also be used to see if the LMM has the same statistical power at the GLMM by setting b1 != 0
b1 <- .2

test<- data.frame(rep=1:nrep)
for(rep in 1:nrep){
  x <- rnorm(nsamples, mean=0, sd=1)
  e <- rnorm(nsamples, mean=0, sd=1)
  beta <- rnorm(ngroups, mean=0, sd=sd.beta)
  x.beta <- rep(beta, each=nsamples.per.group)
  
  # simulation
  Z <- b0 + x * b1 + x.beta + e
  
  p <- inv.logit(Z)
  Y <- rbinom(nsamples, size=1, prob=p)
  
  # These use the summary() function for extracting P-values
  test$p.value.lmer[rep] <- summary(lmer(Y ~ x + (1|group)))$coef[2,5]
  test$p.value.glmer[rep] <- summary(glmer(Y ~ x + (1|group), family=binomial))$coef[2,4]
}

# Histograms of the P-values. Values for the LMM are in gray and for the GLMM in red.
par(mfrow=c(1,1))
hist(test$p.value.lmer, breaks=20, main=paste0("b1 = ",b1,";  lmer: ", mean(test$p.value.lmer < 0.05), " glmer: ", mean(test$p.value.glmer < 0.05)))
hist(test$p.value.glmer, breaks=20, add=T, col="red")

# This shows that the LMM has almost the same power as the GLMM.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#########################################################
# 1.6 Reiteration of results
#########################################################
# Here are all nine models you've run so far

# Route-level methods
summary(lm(asin(sqrt(MEAN_GROUSE)) ~ MEAN_WIND, data=w))$coef

w$SUCCESS <- cbind(w$GROUSE, w$STATIONS - w$GROUSE)
summary(glm(SUCCESS ~ MEAN_WIND, family = binomial, data=w))$coef

summary(glm(SUCCESS ~ MEAN_WIND, family = quasibinomial, data=w))$coef

summary(glmer(SUCCESS ~ MEAN_WIND + (1 | ROUTE), data=w, family=binomial))$coef

# Station-level methods
summary(lm(GROUSE ~ WIND, data=d))$coef

summary(glm(GROUSE ~ WIND, data=d, family=binomial))$coef

summary(glm(GROUSE ~ WIND + ROUTE, data=d, family=binomial))$coef[1:2,]

summary(glmer(GROUSE ~ WIND + (1 | ROUTE), data=d, family=binomial))$coef

summary(lmer(GROUSE ~ WIND + (1 | ROUTE), data=d))$coef


