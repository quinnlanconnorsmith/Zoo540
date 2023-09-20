#Here's code and notes from Tony's book - Correlated Data 
#Errors are independent 
#Error has sme variance 
#yi depends linearly on xi 
#xi and error are uncorrelated 

d <- read.csv(file="grouse_data.csv", header=T)
d
# STATION and ROUTE were uploaded as integers; this converts them to factors.
d$STATION <- as.factor(d$STATION)
d$ROUTE <- as.factor(d$ROUTE)

#Peek at the fist 20 rows of data
head(d, 20)

# Combine all values for each ROUTE.
w <- data.frame(aggregate(cbind(d$GROUSE, d$WIND, d$LAT, d$LONG),
                          by = list(d$ROUTE), FUN = mean))
# For clarity, I've added the column names.
names(w) <- c("ROUTE","MEAN_GROUSE","MEAN_WIND","LAT","LONG")
# Add the count of the number of GROUSE and STATIONS per ROUTE.
w$GROUSE <- aggregate(d$GROUSE, by = list(d$ROUTE), FUN = sum)[,2]
w$STATIONS <- aggregate(array(1, c(nrow(d), 1)),
                        by = list(d$ROUTE), FUN = sum)[,2]
head(w)

#1.4.1
# LM: Analyses at the ROUTE level using an arcsine-square-root transform
summary(lm(asin(sqrt(MEAN_GROUSE)) ~ MEAN_WIND, data=w))

#1.4.2
w$SUCCESS <- cbind(w$GROUSE, w$STATIONS - w$GROUSE)
summary(glm(SUCCESS ~ MEAN_WIND, family = binomial, data=w))

# Likelihood Ratio Test for b1.
mod.f <- glm(SUCCESS ~ MEAN_WIND, family = binomial, data=w)
mod.r <- glm(SUCCESS ~ 1, family = binomial, data=w)
deviance <- 2*(logLik(mod.f) - logLik(mod.r))
pchisq(deviance, df=1, lower.tail=F)
LRT.b1 <- c(dev = deviance, p.value=pchisq(deviance, df=1, lower.tail=F))

#1.4.3
summary(glm(SUCCESS ~ MEAN_WIND, family = quasibinomial, data=w))

#1.4.4
library(lme4)
summary(glmer(SUCCESS ~ MEAN_WIND + (1 | ROUTE), data=w, family=binomial))

#1.4.5
b0 <- 0
b1 <- 0
n <- 8
nsamples <- 1000
x <- rnorm(nsamples, mean=0, sd=1)
# This is set to 0 or 1 for the left and right panels of figure 1.2
sd.e <- 0
e <- rnorm(nsamples, mean=0, sd=sd.e)
Z <- b0 + b1 * x + e
library(boot)
p <- inv.logit(Z)
Y <- rbinom(nsamples, size=n, prob=p)
