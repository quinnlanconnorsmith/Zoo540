#ZOO 540 PS6 
#Cassie, Hangkai, Quinn
# Packages you need
library(lme4)
library(lmerTest)
library(car)
library(tidyverse)
library(nlme)
library(ggiraphExtra)
library(ggiraph)

# This is a function that analyzes a data.frame called "boot" in which each column contains bootstrap values of coefficients (from simulations). It returns the mean, sd, 95% confidence intervals, and P-value for H0:coefficient = 0. It is assumed that the first row of "boot" gives the base parameter estimates from the real data that were used to simulate the data.
boot_coef <- function(base, boot) {
  # construct a data.frame to collect information from "boot"; "boot" is given informative row and column names.
  conf <- data.frame(matrix(NA, nrow = 5, ncol = ncol(boot)))
  names(conf) <- names(boot)
  rownames(conf) <-
    c("mean",
      "se",
      "lower95",
      "upper95",
      "P-value")
  
  for (i.coef in 1:ncol(boot)) {
    vals <- sort(boot[, i.coef])
    conf[1, i.coef] <- mean(vals)
    conf[2, i.coef] <- sd(vals)
    conf[3, i.coef] <- vals[ceiling(.025 * nrow(boot))]
    conf[4, i.coef] <- vals[floor(.975 * nrow(boot))]
    conf[5, i.coef] <- mean(base[i.coef] * vals < 0)
    
  }
  return(t(conf))
}
#################################################################
# Dataset
#################################################################
# There are 2 sites, CO (county line seep) and BS (banana slug seep), and 4 sites total, 3 at CO and 1 at BS. I recorded which plant each seedpod came from (plant_id); there are 82 plants in the study across 399 observations. The response variables are the number of seeds which I have broken up as the number of seeds in the pod (seeds, which is between 0-4), and whether there were seeds present (seed.present, 0 or 1). I think this makes the data sort of cool, because you can look at both the binary yes/no and the amount of seeds as response variables. Then I have the treatments of which there are 8 groups:

# C: control
# U: unbagged
# +: pollen supplementation
# B: bagged
# A: ash

# This is a 2 x 2 x 2 crossed experiment with ash, pollinator exclusion (bag), and pollen supplementation as the variables manipulated. In addition to the treatments category, I have broken up the treatment into its parts in the dataset: 

# ash (ash/control), 
# bag (bag/unbagged)
# pollen (supplement or none). 

# This way you can look at each factor in a model.

# plot.id is coded as 1, 2, 3, so below site.plot.id is created by concatenating site.id and plot.id to give a factor with four levels.

# read data
d <- read.csv(file="Melone seedset data 25Oct23.csv", header=T)

# Do a little data wrangling
d$date <- as.Date(d$date, format = "%m/%d/%y")
d$datecat <- as.factor(d$date)
d$site.id <- as.factor(d$site.id)
d$plant.id <- as.factor(d$plant.id)
d$site.plot.id <- as.factor(paste0(d$site.id, d$plot.id))
d$trt <- as.factor(d$trt)
d$ash <- as.factor(d$ash)
d$bag <- as.factor(d$bag)
d$pollen <- as.factor(d$pollen)

# display data
summary(d)
head(d,20)

# This is just to look at the data by plant.id to see how much data came from each plant
plant.distribution <- aggregate(seed.present ~ plant.id, FUN = length, data = d)
names(plant.distribution)[2] <- "number.of.flowers.per.plant"
plant.distribution

# This is just to look at the data by site.plot.id to see how much data came from each plant
site.plot.distribution <- aggregate(seed.present ~ site.plot.id, FUN = length, data = d)
names(site.plot.distribution)[2] <- "number.of.flowers.per.plant"
site.plot.distribution

# Question: What Gigi wants to know:
# 	Does ash affect seed set? 
# 	How strong is this effect? 
# 	What combination of treatment factors influences seed set? 
# 
# Issues with the data:
# 	
# 1. Hierarchical structure
# 2. Discrete (non-Gaussian) dependent variables (seeds and seed.present)
# 3. Unequal number of plants per site.plot and flowers per plant (although this might not be an issue)
junk <-paste0(d$pollen,d$bag)
d$trt <- paste0(junk,d$ash)

d <- d %>%
  mutate(trt = recode(trt, noneunbaggedcontrol = 'control' ))

#So going forward, I'm treating 'noneunbaggedcontrol' as the true control group

ggplot(data=d, aes(x = ash, y = seeds)) +
  geom_boxplot(outlier.size=2, outlier.shape=21)+
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Ash seems to have some influence on seed count in ashed plants 

ggplot(data=d, aes(x = trt, y = seeds, color=trt)) +
  geom_boxplot(outlier.size=2, outlier.shape=21)+
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Model selection process

model1 <- lm(seeds~trt, data=d)
summary(model1)
plot(model1, ask=F)

model2 <- lm(seeds~trt + site.id, data=d)
summary(model2)
plot(model2, ask=F)

anova(model1,model2)

#Site has a slightly significant influence on # of seeds

model3<- lm(seeds~trt + plant.id, data=d)

anova(model2,model3)


model4 <- lm(seeds~trt + site.id + plant.id, data=d)
anova(model2,model4) 
anova(model3,model4)  
#Plant ID is contributing most of the sig

model5 <- lm(seeds~trt + plant.id + trt*plant.id, data=d)

anova(model3,model5)  
#Interaction is not significant 

#So lets stick with model 3
model3 <- lm(seeds~trt + plant.id, data=d)
#plot(d$site.id, residuals(model3, type='pearson'))
plot(d$plant.id, residuals(model3, type='pearson'))
abline(h=0)
#Moving this to an lme random effect and double checking

model6 <- lme(seeds~trt, random= ~1|site.id, method = 'REML', data=d)
summary(model6)
plot(model6)


model7 <- lme(seeds~trt, random= ~1|plant.id, method = 'REML', data=d)
summary(model7)
plot(model7)

model8 <- lme(seeds~trt, random= ~1|plant.id/site.id, method = 'REML', data=d)
summary(model8)
plot(model8)

anova(model6,model7)
anova(model7,model8)

#model 7 - random effect of plant

#Refit model with ML
#Don't need to do this, but would need to if there was an interaction we wanted to test

#model9 <- lme(seeds~trt, random= ~1|plant.id, method = 'ML', data=d)

seed_model_final <-lme(seeds~trt, random= ~1|plant.id, method = 'REML', data=d)
hist(residuals(seed_model_final))
seed_final_residuals <-residuals(seed_model_final, type='pearson')
seed_final_fitted <-fitted.values(seed_model_final)
#Checking assumptions 
plot(seed_final_fitted,seed_final_residuals)
abline(h=0)

#plot(seed_final_residuals~trt, data=d)
#abline(h=0)
plot(seed_final_residuals~plant.id, data=d)
abline(h=0)

summary(seed_model_final)

var_rep_resid <- VarCorr(seed_model_final)
var_rep <- as.numeric(var_rep_resid[1])
var_resid <- as.numeric(var_rep_resid[2])

var_rep/(var_rep + var_resid)
#Correlation is low between obs (0.19)

#Visualizing the GLMM with random effect of plant.id

sjPlot::plot_model(seed_model_final)

sjPlot::plot_model(seed_model_final, 
                   show.values=TRUE, show.p=TRUE,
                   title="Effect of Treatments Seed Count")

sjPlot::tab_model(seed_model_final, 
                  show.re.var= TRUE, 
                  dv.labels= "Effect of Treatments on Seed Count")

seed_pred <- expand.grid(
  seeds = seq(0,4,1),
  trt = d$trt, 
  plant.id = levels(d$plant.id)
)
seed_pred

seed_pred$seeds <-predict(seed_model_final, newdata = seed_pred)

#Predictions for each treatment from all the plants - these are not values from the sjplot
#might not be as useful because it's tough to predict data from plants that did not receive all the treatments
  
ggplot(data=seed_pred, aes(x = trt, y = seeds)) +
  geom_boxplot(outlier.size=2, outlier.shape=21)+
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
####Goofin####
#Just looking at effect of Ash here keeping plant.id as random. If you actually want to do this go back through the model selection! 

seed_ash_model <-lme(seeds~ash, random= ~1|plant.id, method = 'REML', data=d)
hist(residuals(seed_ash_model))
seed_ash_residuals <-residuals(seed_ash_model, type='pearson')
seed_ash_fitted <-fitted.values(seed_ash_model)
#Checking assumptions 
plot(seed_ash_fitted,seed_ash_residuals)
abline(h=0)

#plot(seed_final_residuals~trt, data=d)
#abline(h=0)
plot(seed_ash_residuals~plant.id, data=d)
abline(h=0)

summary(seed_ash_model)

var_rep_resid <- VarCorr(seed_model_final)
var_rep <- as.numeric(var_rep_resid[1])
var_resid <- as.numeric(var_rep_resid[2])

var_rep/(var_rep + var_resid)

sjPlot::plot_model(seed_ash_model, 
                   show.values=TRUE, show.p=TRUE,
                   title="Effect of Treatments Seed Count")

