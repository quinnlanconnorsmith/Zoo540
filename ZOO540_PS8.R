#Cassie, Hangkai, Quinn
# Packages you might need
library(lme4)
library(lmerTest)
library(car)

#################################################################
# Dataset
#################################################################

# To see all of the text below, remember to use the "wrap lines" option in the Preferences.

# DATA
# "zoo540_images_2023.csv" contains the data relevant to landscape preference. Each respondent (n = 402) viewed & stated their preference for four images - one respondent (R_1EnpXMq6Eq1QHi8) skipped their first image, thus one row is missing (total rows 402 x 4 - 1 = 1607). 

# ResponseId is an identifier for each survey participant 
# Image_Set is the image set used in the survey (1, 2, or 3)
# Image is the image within the image set (1, 2, 3, or 4)
# Location identifies where the photo was taken. Each location has an RCP 4.5 & an RCP 8.5 edit
# Preference is the preference for one scene over the other (-3 to 3), with -3 indicating strong preference for projection and 3 indicating strong preference for current scene (2022)
# RCP is the climate scenario the edited image represents
# Scenic is whether the image contains additional scenic elements (e.g. mountain vista, geyser)

# MR1_Experience are factor scores relating to experiences (Exploratory factor analysis with varimax rotation).
# MR2_Access are factor scores relating to access
# MR3_Engagement are factor scores relating to engagement

# read raw data
d.raw <- read.csv(file="Zoo540_PS8_images_2023.csv", header=T)
summary(d.raw)

# two respondents skipped the section relevant for factor scores. these responses could be excluded. 
d <- d.raw[d.raw$ResponseId != "R_814YmULtWKDPVZL" & d.raw$ResponseId != "R_jzqdqU1bYiXy42w",]

# More or less useful transformations for 1) scenic column, and 2) creating preference groups
d$Scenic <- ifelse(d$Scenic == "N", 0, 1)
d$Preference.Simple <- ifelse(d$Preference < 0, "Future", 
                              ifelse(d$Preference > 0, "Present", "None"))

# convert to factors
d$ResponseId <- as.factor(d$ResponseId)
d$Image.Set <- as.factor(d$Image.Set)
d$Image <- as.factor(d$Image)
d$Location <- as.factor(d$Location)
d$RCP <- as.factor(d$RCP)
d$Scenic <- as.factor(d$Scenic)
d$Preference.Simple <- as.factor(d$Preference.Simple)

summary(d)
# Questions: 
# 1) Do RCP and Scenic Elements affect visitor preference for future landscapes? 
# 2) Do factor scores affect visitor preference? 

########################################################################
# visualizations
########################################################################
aggregate(ResponseId ~ Preference, data = d, FUN = length)
hist(d$Preference-.01)

MR1 <- log()

library(tidyverse)

#Just checking some visuals 
par(mfrow = c(1,1))
hist((d$MR1.Experience))
hist((d$MR2.Access))
hist((d$MR3.Engagement))

#hist(sqrt(d$MR1.Experience))
#hist(sqrt(d$MR2.Access))
#hist(sqrt(d$MR3.Engagement))

plot(d$MR1.Experience~d$Preference)
abline(lm(d$MR1.Experience~d$Preference))
plot(d$MR2.Access~d$Preference)
abline(lm(d$MR2.Access~d$Preference))
plot(d$MR3.Engagement~d$Preference)
abline(lm(d$MR3.Engagement~d$Preference))

#all lm seem to have similar slopes and intercepts 

ggplot(d, aes(x=Preference)) +
  geom_histogram() +
  facet_wrap(~RCP)

ggplot(d, aes(x=Preference)) +
  geom_histogram() +
  facet_wrap(~Scenic)

ggplot(d, aes(x=Preference)) +
  geom_histogram() +
  facet_grid(RCP~Scenic)

#So lots of folks said 3 no matter the RCP or whether there were scenic elements 
#Let's feed this all into a model and see what we can find 

#Starting off with a simple lm 

mod1 <- lm(Preference ~ RCP + Scenic, data=d)
summary(mod1)
mod2 <- lm(Preference ~ RCP + Scenic + MR1.Experience + MR2.Access + MR3.Engagement, data=d)
summary(mod2)

anova(mod1,mod2)

#So in the linear models, the M1-M3 variables are important for explaining preference 

#Because there's multiple data points for people, let's add a random effect for responseid
ggplot(d, aes(x=Preference.Simple)) +
  geom_histogram(stat="count")

#mod3 <- glmer(Preference.Simple ~ RCP + Scenic + MR1.Experience + MR2.Access + MR3.Engagement + (1|ResponseId), family=binomial, data=d)

#summary(mod3)

#NOT doing this, as preference.simple isn't binomial

#So here, we lose information going from the preference to binomial preference.simple
#BUT scenic comes out as significant? 

#Let's go with a lmer 
mod4 <- lmer(Preference ~ RCP + Scenic + MR1.Experience + MR2.Access + MR3.Engagement + (1|ResponseId), data=d)

summary(mod4)

#What about interactions? 

mod5 <- lmer(Preference ~ RCP + Scenic + MR1.Experience + MR2.Access + MR3.Engagement + RCP*Scenic  + (1|ResponseId), data=d)
summary(mod5)

anova(mod4,mod5)
anova(mod5)

#So there is an interaciton between RCP scenario and Scenic choices
#How does this work with 4.5 vs. 8.5? are we losing the 4.5?
#I think this model is "telling" us that in the "worse" climate scenario, scenic images are different 

mod6 <- lmer(Preference ~ RCP + Scenic + MR1.Experience + MR2.Access + MR3.Engagement + RCP*Scenic  + MR1.Experience*MR2.Access*MR3.Engagement +(1|ResponseId), data=d)
summary(mod6)

#Here's an overfit model that also has an interaction between MR1-MR3. It seems that the interaction between Access and Engagement may be significant

# Another approach, more models 
## simple linear model - RCP important
lm1<-lm(Preference ~ RCP, data=d)
summary(lm1)
par(mfrow=c(2,2))
plot(lm1)
### scenic not as important
lm2<-lm(Preference ~ Scenic, data=d)
summary(lm2)
par(mfrow=c(2,2))
plot(lm2)

# LMER1
# mixed model with RCP and Scenic as fixed and response id as random
# response id .4672 variance
# AIC 4798.2
lmer1<-lmer(Preference ~ RCP + Scenic  + (1|ResponseId), REML=FALSE, data=d)
summary(lmer1)
par(mfrow=c(1,1))
plot(lmer1)

# LMER2
# response id .3468 lower than before when experience is added
# AIC 4793.5
lmer2<-lmer(Preference ~ RCP + (1 + MR1.Experience |ResponseId), REML=FALSE, data=d)
summary(lmer2)
par(mfrow=c(1,1))
plot(lmer2)


# engagement does not improve model
lmer3<-lmer(Preference ~ RCP + (1 + MR3.Engagement |ResponseId), REML=FALSE, data=d)
summary(lmer3)
par(mfrow=c(1,1))
plot(lmer3)

# LMER4
# AIC 4789.2
# response id .3481 experience .1632
# interaction of RCP and scenic improves model
lmer4<-lmer(Preference ~ RCP * Scenic + (1 + MR1.Experience |ResponseId), REML=FALSE, data=d)
summary(lmer4)
par(mfrow=c(1,1))
plot(lmer4)

anova(lmer3, lmer4, test="Chisq")

# LMER5
lmer5 <- lmer(Preference ~  RCP + Scenic + MR1.Experience*MR2.Access +
                MR1.Experience*MR3.Engagement + (1|ResponseId), data = d, REML = FALSE)
summary(lmer3)

anova(lmer1, lmer3, test="Chisq")

#Round 3 - hybrid 

# Transform 'Scenic' into a binary variable and create 'Preference.Simple' for analysis
d$Scenic <- ifelse(d$Scenic == "N", 0, 1)
d$Preference.Simple <- ifelse(d$Preference < 0, "Future", ifelse(d$Preference > 0, "Present", "None"))

# Convert columns to factors for analysis
d$ResponseId <- as.factor(d$ResponseId)
d$Image.Set <- as.factor(d$Image.Set)
d$Image <- as.factor(d$Image)
d$Location <- as.factor(d$Location)
d$RCP <- as.factor(d$RCP)
d$Scenic <- as.factor(d$Scenic)

# Exploratory Data Analysis
# Histograms and aggregate plots to understand data distribution
hist(d$Preference)
aggregate(ResponseId ~ Preference, data = d, FUN = length)

# Question 1 Analysis: Do RCP and Scenic Elements affect visitor preference for future landscapes?
# Mixed model with RCP and Scenic as fixed effects and ResponseId as a random effect
lmer1 <- lmer(Preference ~ RCP + Scenic + (1|ResponseId), REML=FALSE, data=d)
summary(lmer1)

# Question 2 Analysis: Do factor scores affect visitor preference?
# Linear mixed model incorporating MR1, MR2, and MR3 scores as fixed effects
# ResponseId as a random effect
mod4 <- lmer(Preference ~ MR1.Experience + MR2.Access + MR3.Engagement + (1|ResponseId), data=d)
summary(mod4)

# Visualizations to check model diagnostics
par(mfrow=c(2,2))
plot(lmer1)
plot(mod4)

