#Quinn Connor & Janelle 
# PS3
# due 3 October, 2023 by midnight

# This problem set is to check to make sure you are up and running in R, and it begins to explore the ruffed grouse data. It is based on correlated_data Ch1 ("Ch_1_10Aug18.R") sections 1.5 to 1.7.

# For the homework, EACH INDIVIDUAL SHOULD TURN IN A SEPARATE .R FILE. Start with the code from correlated_data ("Ch_1_10Aug18.R") and add to it anything you need. Identify you new code (so I can find it) by placing it between marker rows #~~~~~~~~~~~~~~~~~~~~~~~~~~~. I will get your full answers to the questions by randomly asking people in class, so there is no need for you to write them down.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. What questions do you have about the material in section 1.5 to 1.7? What needs more explanation? I'm serious about asking this question, because I want to improve the book. (NOTE: This requires NO NEW R CODE.)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Simulations seem like a great way to explore modeling methodologies and take a deeper look at what we're doing. What situations does simulating not make sense, if any? Is it a waste of time in these situations, or detrimental (e.g. leading down the "wrong" path)? 

#What is the difference, in your words, between anova() and Anova()? 
?anova
?Anova
#Ope you cover this in correlated data - disregard 
#For the grouse data, 50 means must be estimated. The analysis is easy to do with the glm() function, although here I have also used the function Anova() from the car library (Fox and Weisberg 2011) because it gives a summary statistic for the overall effect of the factor ROUTE in the model, rather than reporting all 49 estimates of b2[factor] for the differences among routes. (The remaining route serves as the intercept, accounting for the total of 50 levels of the factor ROUTE.)

#Would love to go over simulation of the lmer with the class! 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Suppose Y1, Y2, and Y3 are independent normal random variables with variances 1, 2, and 3, respectively. You can think of (Y1 + Y3) and (Y2 + Y3) as the number of grouse observed in station 1 and station 2 of route 1: Y3 gives the random effect that is common to all stations within route 1, while Y1 and Y2 give the station-specific random effects for stations 1 and 2 within route 1. What is the covariance between (Y1 + Y3) and (Y2 + Y3)?  If you want, check your answer with a simulation.

#  Cov(Y1 + Y3, Y2 + Y3) = Cov(Y1, Y2) + Cov(Y1, Y3) + Cov(Y2, Y3) + Cov(Y3, Y3) Since Y1, Y2, and Y3 are independent, Cov(Y1, Y2) = Cov(Y1, Y3) = Cov(Y2, Y3) = 0. Also, the covariance of a random variable with itself is equal to its variance, so Cov(Y3, Y3) = Var(Y3) = 3.The covariance between (Y1 + Y3) and (Y2 + Y3) simplifies to: Cov(Y1 + Y3, Y2 + Y3) = 0 + 0 + 0 + 3 = 3                    The covariance between (Y1 + Y3) and (Y2 + Y3) is 3
#From Tony on Monday
predicted.cov <-3
Y1=rnorm(1000, mean=0, sd=1)
Y2=rnorm(1000, mean=0, sd=sqrt(2))
Y3=rnorm(1000, mean=0, sd=sqrt(3))
Y4=Y1+Y3
Y5 =Y2+Y3
cov(Y4,Y5)
predicted.cov
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Anova(glm(GROUSE ~ WIND + ROUTE, data=d.no0, family=binomial))
summary(glmer(GROUSE ~ WIND + (1 | ROUTE), data=d.no0, family=binomial))

######
# 1.5.5 LMM: Even though this is a linear model, it incorporates the non-independence of STATIONS within ROUTES.
summary(lmer(GROUSE ~ WIND + (1 | ROUTE), data=d))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3.  (from book) Fit a LM to the station-level data treating route as a factor. This is similar to the LMM treating route as a random effect (subsection 1.5.5). 1). Do they give similar results about the effect of WIND on GROUSE? 2).Is the comparison between these LM and LMM similar to the comparison between the GLM treating route as a factor (subsection 1.5.3) and the GLMM treating route as a random effect (subsection 1.5.4)? 3). Produce a figure like figure 1.3. Does this explain your answer in the same way figure 1.3 explains the difference between the GLM and GLMM?

#1). Where the wind is very significant in the LMM, it's only close to significant in the LM
summary(lm(GROUSE ~ WIND + ROUTE, data=d))
summary(lmer(GROUSE ~ WIND + (1 | ROUTE), data=d))

lm <-lm(GROUSE ~ WIND + ROUTE, data=d)
lmer <-lmer(GROUSE ~ WIND + (1 | ROUTE), data=d)
#2). From above
z.glm <- glm(GROUSE ~ 0 + WIND + ROUTE, data=d, family=binomial)
z.glmm <- glmer(GROUSE ~ WIND + (1 | ROUTE), data=d, family=binomial)
summary(z.glm)
summary(z.glmm)
#In a way, yes. The glmm which treats route as a random effect has a significant effect of wind, where the glm does not. Do the glm and lm "spread" out our power for detecting statistical significance? 

#3). We don't have that far grouping like we saw with the glm and glmm, instead there's just a slight spread in the LM. Only positive values in the lmer? 
# Fig. 1.3
par(mfrow=c(1,2), mai=c(1,.8,.5,.1), mgp=c(3,.5,0))

# This plots the 50 values of route-factor coefficients
hist(lm$coef[2:51], breaks=(-5:5), freq=F,  main="Estimates from LM", xlab="b1")

# This plots the 50 values of the random effects associated with the 50 routes
hist(coef(lmer)$ROUTE[[1]], breaks=(-5:5), freq=F, main="Estimates from lmer", xlab="b1")
lines(.02*(-1000:1000), dnorm(.02*(-1000:1000), mean=fixef(lmer)['WIND'], sd=summary(lmer)$var[[1]][1]^.5), col="black")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. In section 1.4.5 there is a simulation for the GLMM model of route-level data. How could you use a model (it has to be a different model) for station-level data to test whether the P-values given by the LMM (section 1.5.5) are correct? This question is similar to Question 4 in Problem Set 2.

#Attempt 1

#sub_n  <- 50 # number of subjects in this simulation
#sub_sd <- 100 # SD for the subjects' random intercept

#sub <- tibble(
#  sub_id = 1:sub_n,
#  sub_i  = rnorm(sub_n, 0, sub_sd), # random intercept
#  sub_cond = rep(c("wind","no_wind"), each = sub_n/2) # between-subjects factor
#)

#Attempt 2

nroute = 50
nstation = 8
mu = 5
sds = 2
sd = 1

route = rep(LETTERS[1:nroute], each = nstation) 
station = rep(LETTERS[1:nstation], each = nroute) 
stationeff = rnorm(nstation, 0, sds) 
stationeff = rep(standeff, each = nroute) 
routeeff = rnorm(nstation*nroute, 0, sd) 
dat = data.frame(station, standeff, route, routeeff) 
dat$resp = with(dat, mu + stationeff + routeeff ) 

#Unsure what to do about this

# Strategy:
# 1. Set up a simulatmodel that can produce data from the same statistical model you are going to test
# 2. Simulate lots of datasets (2000)
# 3. For each data set, fit the model.
# 4. Compare the results of the large numbers of fits to the model that you used to simulate the data
#summary(lm(GROUSE ~ WIND + ROUTE, data=d))
#lm(GROUSE ~ WIND + ROUTE, data=d)$coef[2,4]
#summary(lmer(GROUSE ~ WIND + (1 | ROUTE), data=d))

# Set up test parameter values that you will estimate with the statistical model 
#b0 <- 0
#b1 <- 0
#n <- 8
#sd.e <- 1

# Set the size of each dataset
#nsamples <- 10

# This just gives an example simulation that you should check.
#x <- rnorm(nsamples, mean=0, sd=1)
#e <- rnorm(nsamples, mean=0, sd=sd.e)

#Z <- b0 + b1 * x + e
#p <- inv.logit(Z)
#Y <- rbinom(nsamples, size=n, prob=p)

summary(lmer(GROUSE ~ WIND + (1 | ROUTE), data=d))

# Perform many simulations to create lots of datasets and fit the model to each.
#nrep <- 2000
#test <- data.frame(rep=1:nrep)
#for(rep in 1:nrep){
#  x <- rnorm(nsamples, mean=0, sd=1)
#  e <- rnorm(nsamples, mean=0, sd=sd.e)
  
#  Z <- b0 + b1 * x + e
#  p <- inv.logit(Z)
#  Y <- rbinom(nsamples, size=n, prob=p)
  
#  test$p.value[rep] <- summary(lmer(GROUSE ~ WIND + (1 | ROUTE), data=d))$coef[2,4]
#  test$b1[rep] <- summary(lmer(GROUSE ~ WIND + (1 | ROUTE), data=d))$coef[2,1]
#}

# Plot a histogram of the P-values and compute the number of datasets for which P < 0.05
#par(mfrow = c(1,1), mai = c(1,.7,1,.3))
#hist(test$p.value, breaks=20, main=paste0("Proportion of P-values < 0.05 = ",mean(test$p.value < 0.05)))

# Plot a histogram of the estimates of b1
#hist(test$b1, breaks=20, main="b1")

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