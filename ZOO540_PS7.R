#Emily, Frank, Quinn
# Packages you might need
library(lme4)
library(lmerTest)
library(car)

#################################################################
# Dataset
#################################################################

# There are two files. To see all of the text below, remember to use the "wrap lines" option in the Preferences.

# RAW DATA
# "qs_zoo540_raw_data_8Nov23.csv" contains all the digitized data from the trip, with 1 row for each measured burrow (991 total)
# date is the date (10/16/2022-10/21/2022) 
# site is the set of treatment plots (1-10) each site had 1 of each treatment
# trt is the treatment that plot received (either Control, Compact, Marine.Debris, or Rake) 
# burrow.diameter is that individual burrow diameter measured at the widest point (a NA is used for plots with no burrows)
# mean.temp is the average air temperature for the island from that day (taken from UGAMI-LTER stations)

# SUMMARY DATA
# "qs_zoo540_sum_data_8Nov23.csv" is the summarized version of the raw data
# The only new things are count (number of burrows in that plot) and mean.burrow.diameter (AVERAGE burrow diameter for that plot) 


# read raw data
d.raw <- read.csv(file="qs_zoo540_raw_data_8Nov23.csv", header=T)
d.raw$date <- as.Date(d.raw$date, format = "%m/%d/%y")
d.raw$trt <- as.factor(d.raw$trt)
summary(d.raw)

# construct a data.frame aggregated to the date/site/plot with mean.burrow.diameter and count
d.agg <- aggregate(cbind(burrow.diameter, mean.temp) ~ date + site + trt, data = d.raw, FUN = mean, drop = F)
d.agg <- cbind(d.agg, aggregate(burrow.diameter ~ date + site + trt, data = d.raw, FUN = length, drop = F)[,4])
names(d.agg)[c(4, 6)] <- c("mean.burrow.diameter", "count")
d.agg$count[is.na(d.agg$count)] <- 0
d.agg$mean.burrow.diameter[is.na(d.agg$mean.burrow.diameter)] <- 0
summary(d.agg)

########################################
# check against "qs_zoo540_sum_data_8Nov23.csv"

# (It is ALWAYS a good idea to meticulously check data before you start. Quinn provided two files, so I decided to check them against each other.)

d <- read.csv(file="qs_zoo540_sum_data_8Nov23.csv", header=T)
d$date <- as.Date(d$date, format = "%m/%d/%y")
d$trt <- as.factor(d$trt)
summary(d)

# merge aggregated fine from d.raw and d
d.merged <- merge(d, d.agg, by = c("date", "site", "trt"))

# These two checks found that there were two rows in the d and d.agg that differ, with d.agg (i.e., the raw data) containing one additional 
d.merged[d.merged$count.x != d.merged$count.y,]
d.merged[abs(d.merged$mean.burrow.diameter.x - d.merged$mean.burrow.diameter.y) > .001 & !is.na(d.merged$mean.burrow.diameter.x),]

# The raw data look okay.
d.raw[d.raw$date == "2022-10-20" & d.raw$site == 2 & d.raw$trt == "Control",]
d.raw[d.raw$date == "2022-10-21" & d.raw$site == 5 & d.raw$trt == "Compact",]

# I checked with Quinn, and he said use the dataset generated from d.raw
d <- d.agg

library(ggplot2)
library(tidyverse)
library(lubridate)
library(RColorBrewer)
library(wesanderson)
library(nlme)
library(leaflet)
library(lme4)
library(ggiraph)
library(ggiraphExtra)
library(ggeffects)
library(lmeSplines)
library(sjPlot)
library(effects)
library(merTools)
library(AICcmodavg)
library(zoo)


d.avg <- aggregate(cbind(count, mean.temp, mean.burrow.diameter) ~ date + trt, data = d.agg, FUN = mean, drop = F)


#Just some quick visuals before we get into it 
ggplot(data=d.agg, aes(x = trt, y = count)) +
  facet_wrap(~date) +
  geom_boxplot(outlier.size=2, outlier.shape=21)+
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(d.agg[which(d.agg$mean.burrow.diameter>0),], aes(x = trt, y = mean.burrow.diameter)) +
  facet_wrap(~date) +
  geom_boxplot(outlier.size=2, outlier.shape=21)+
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Adding a second axis (TEMP)
#Burrow count 
ggplot() +
  geom_point(data=d.avg, aes(x = date, y = count, color=trt), size=5) +
  geom_line(data=d.avg, aes(x = date, y = count, color=trt, group=trt), linewidth=1) +
  geom_point(data=d.agg, aes(x = date, y = mean.temp/1.5), size=2) +
  geom_line(data=d.agg, aes(x = date, y = mean.temp/1.5, color="Temperature", group=trt), linewidth=1) +
  geom_jitter(data=d.agg, aes(x = date, y = count, color=trt), width=0.1, alpha=0.5) +
  scale_y_continuous(name="Number of Burrows", sec.axis = sec_axis(~.*1.5, name="Average Air Temperature (°C)")) +
  labs(x="Date") +
  guides(color=guide_legend(title="Treatment"))+
  scale_colour_manual(values = c("#56B4E9", "#999999", "#009E73", "#D55E00","black")) +
  theme(text = element_text(size = 13)) 

#Burrow diameter
ggplot() +
  geom_point(data=d.avg, aes(x = date, y = mean.burrow.diameter, color=trt), size=5) +
  geom_line(data=d.avg, aes(x = date, y = mean.burrow.diameter, color=trt, group=trt), linewidth=1) +
  geom_point(data=d.agg, aes(x = date, y = mean.temp*2.5), size=2) +
  geom_line(data=d.agg, aes(x = date, y = mean.temp*2.5, color="Temperature", group=trt), linewidth=1) +
  geom_jitter(data=d.agg, aes(x = date, y = mean.burrow.diameter, color=trt), width=0.1, alpha=0.5) +
  scale_y_continuous(name="Mean Burrow Diameter (mm)", sec.axis = sec_axis(~./2.5, name="Average Air Temperature (°C)"))+
  labs(x="Date")+
  guides(color=guide_legend(title="Treatment"))+
  scale_colour_manual(values = c("#56B4E9", "#999999", "#009E73", "#D55E00","black")) +
  theme(text = element_text(size = 13)) 

#Converting data to comp data 

d.comp <- d.agg %>%
  # Convert the date column to a Date object if it's not already
  mutate(date = as.Date(date)) %>%
  # Arrange the data by date, site, and treatment
  arrange(date, site, trt) %>%
  # Group by site and treatment
  group_by(site, trt) %>%
  # Calculate the difference in counts between each date and the first date
  mutate(count.diff = count - first(count)) %>%
  mutate(dia.diff = mean.burrow.diameter-first(mean.burrow.diameter))


d.comp2 <- d.comp %>%
  arrange(date) %>%
  filter(row_number() > 1)

write.csv(d.comp2, "d_comp2.csv", row.names=FALSE)
#I need to change the mean.temp NA's to the value from that day and I can't do that in R, so I wrote the .csv and did it manually 

d.comp3 <- d_comp2 %>%
  # Use mutate to replace a specific value in the column
  mutate(trt = ifelse(trt == "Control", "aControl", trt)) %>%
  mutate(trt = ifelse(trt == "Marine.Debris", "Marine_debris", trt))


#d.comp3 <- d.comp2 %>%
#  mutate(date = as.Date(date)) %>%
  # Arrange the data by date
#  arrange(date, site) %>%
  # Use na.locf to fill missing values with the last observed value
#  fill(mean.temp)

#Modeling

cmod1 <- lm(count.diff~site, data=d.comp3)
summary(cmod1)
plot(cmod1, ask=F)
#Variance of burrow count comp appears homogeneous among replicates 

cmod2 <- lm(count.diff~site+mean.temp, data=d.comp3)
anova(cmod1, cmod2)

cmod3 <-lm(count.diff~trt + mean.temp + trt*mean.temp, data=d.comp3)
summary(cmod3)
plot(d.comp2$site, residuals(cmod3, type='pearson'))

cmod4 <- lme(count.diff~trt + mean.temp + trt*mean.temp, random= ~1|site, method = 'REML', data=d.comp3)
summary(cmod4)
plot(cmod4)

gls_cmod1 <- gls(count.diff~trt + mean.temp + trt*mean.temp, data=d.comp3)
anova(gls_cmod1, cmod4)

#So the random effect can be kept! 
#Checking standardized residials and fitted values 

cmod4_residuals <- residuals(cmod4, type='pearson')
cmod4_fitted <- fitted.values(cmod4)
plot(cmod4_fitted,cmod4_residuals)
abline(h=0)
#plot(cmod4_residuals~trt, data=d.comp3) Unsure why trt makes it angry 
abline(h=0)
plot(cmod4_residuals~mean.temp, data=d.comp3)
abline(h=0)
plot(cmod4_residuals~site, data=d.comp3)
abline(h=0)

#Refit model with ML and drop interaction 
cmod5 <- lme(count.diff~trt + mean.temp + trt*mean.temp, random= ~1|site, method = 'ML', data=d.comp3)
cmod6 <- lme(count.diff~trt + mean.temp, random= ~1|site, method = 'ML', data=d.comp3)
anova(cmod5,cmod6)
#So the interaction term can be dropped 

cmod7 <- lme(count.diff~trt + mean.temp, random= ~1|site, method = 'ML', data=d.comp3)
cmod7_droptreat <- lme(count.diff~mean.temp, random= ~1|site, method = 'ML', data=d.comp3)
cmod7_droptemp <- lme(count.diff~trt, random= ~1|site, method = 'ML', data=d.comp3)

summary(cmod7)
summary(cmod7_droptreat)
summary(cmod7_droptemp)

anova(cmod7,cmod7_droptemp)
anova(cmod7,cmod7_droptreat)
#So both of these are very significant 

summary(cmod7)

par(mfrow = c(1,1))

count.mod.final <-lme(count.diff~trt + mean.temp, random= ~1|site, method = 'REML', data=d.comp3)
hist(residuals(count.mod.final), main="A")
count.mod.residuals <-residuals(count.mod.final, type='pearson')
count.mod.fitted <-fitted.values(count.mod.final)
#Checking assumptions 
qqnorm(count.mod.final, main="B")
#qqline(count.mod.final)

plot(count.mod.fitted,count.mod.residuals, main="C")
abline(h=0)

#The rest 
#plot(count.mod.residuals~trt, data=d.comp3, main="D") Again unsure why trt makes it angry
abline(h=0)
plot(count.mod.residuals~mean.temp, data=d.comp3, main="E")
abline(h=0)
plot(count.mod.residuals~site, data=d.comp3, main="F")
abline(h=0)

summary(count.mod.final)

count.var.resid <- VarCorr(count.mod.final)
a <- as.numeric(count.var.resid[1])
b <- as.numeric(count.var.resid[2])

a/(a + b)

#Correlation between observations is low (0.24)


#Test to see if random slopes were needed 
count.mod.final.slope <-lme(count.diff~trt + mean.temp, random= (~1+mean.temp|site), method = 'REML', data=d.comp3)

anova(count.mod.final, count.mod.final.slope)

#Random slope is not needed

#Burrow diameter time with new comp response 

cmod8 <- lm(dia.diff~site, data=d.comp3)
summary(cmod8)
plot(cmod8, ask=F)
#Varaince of burrow diameter comp appears homogeneous among replicates 

cmod9 <- lm(dia.diff~site+mean.temp, data=d.comp3)
anova(cmod8, cmod9)

cmod10 <-lm(dia.diff~trt + mean.temp + trt*mean.temp, data=d.comp3)
summary(cmod10)
plot(d.comp3$site, residuals(cmod3, type='pearson'))

cmod11 <- lme(dia.diff~trt + mean.temp + trt*mean.temp, random= ~1|site, na.action=na.omit, method = 'REML', data=d.comp3)
summary(cmod11)
plot(cmod11)

gls_cmod2 <- gls(dia.diff~trt + mean.temp + trt*mean.temp, na.action = na.omit, data=d.comp3)
anova(gls_cmod2, cmod11)

#So the random effect can be kept! But it's real close 
#Checking standardized residuals and fitted values 

cmod11_residuals <- residuals(cmod11, type='pearson')
cmod11_fitted <- fitted.values(cmod11)
plot(cmod11_fitted,cmod11_residuals)
abline(h=0)
#plot(cmod11_residuals~trt, data=d.comp3) Angry because of NA's
abline(h=0)
plot(cmod11_residuals~mean.temp, data=d.comp3)
abline(h=0)
plot(cmod11_residuals~site, data=d.comp3)
abline(h=0)

#Refit model with ML and drop interaction 
cmod12 <- lme(dia.diff~trt + mean.temp + trt*mean.temp, random= ~1|site, method = 'ML', data=d.comp3)
cmod13 <- lme(dia.diff~trt + mean.temp, random= ~1|site, method = 'ML', data=d.comp3)
anova(cmod12,cmod13)
#So the interaction term can be dropped - but just barely 

cmod14 <- lme(dia.diff~trt + mean.temp, random= ~1|site, method = 'ML', data=d.comp3)
cmod14_droptreat <- lme(dia.diff~mean.temp, random= ~1|site, method = 'ML', data=d.comp3)
cmod14_droptemp <- lme(dia.diff~trt, random= ~1|site, method = 'ML', data=d.comp3)

summary(cmod14)
summary(cmod14_droptreat)
summary(cmod14_droptemp)

anova(cmod14,cmod14_droptemp)
anova(cmod14,cmod14_droptreat)
#So both of these are significant - treatment is barely so 

summary(cmod14)

dia.model.final <-lme(dia.diff~trt + mean.temp, random= ~1|site, method = 'REML', data=d.comp3)
hist(residuals(dia.model.final), main="A")
dia.model.residuals <-residuals(dia.model.final, type='pearson')
dia.model.fitted <-fitted.values(dia.model.final)
#Checking assumptions 
qqnorm(dia.model.residuals, main="B")
#qqline(dia.model.final)
plot(dia.model.fitted,dia.model.residuals, main="C")
abline(h=0)

#The rest 
#plot(dia.model.residuals~trt, data=d.comp3, main="D")
abline(h=0)
plot(dia.model.residuals~mean.temp, data=d.comp3, main="E")
abline(h=0)
plot(dia.model.residuals~site, data=d.comp3, main="F")
abline(h=0)

summary(dia.model.final)

ctrl <- lmeControl(opt='optim')
dia.model.final.slope <-lme(dia.diff~trt + mean.temp, random= ~1+trt|site, control=ctrl, method = 'REML', data=d.comp3)

anova(dia.model.final.slope, dia.model.final)

summary(dia.model.final)
#Compact is sig dif 
#Rake is not sig diff 
#MD is not sig dif 
#Treamtent is actually significant (barely) when looking at comparison 

dia.rep.resid <- VarCorr(dia.model.final)
c <- as.numeric(dia.rep.resid[1])
d <- as.numeric(dia.rep.resid[2])

c/(c+d)

#Correlation between observations is low (0.04)

sjPlot::plot_model(count.mod.final)

sjPlot::plot_model(count.mod.final,
                   show.values=TRUE, show.p=TRUE)
sjPlot::plot_model(count.mod.final, 
                   axis.labels=c("Temperature", "Rake", "Marine Debris", "Compact"),
                   show.values=TRUE, show.p=TRUE,
                   title="Effect of Treatments on Burrow Count")

sjPlot::tab_model(count.mod.final, 
                  show.re.var= TRUE, 
                  pred.labels =c("(Intercept)","Compact", "Marine Debris", "Rake", "Temperature"),
                  dv.labels= "Treatment and Temperature Effects on Burrow Density")

sjPlot::tab_model(dia.model.final, 
                  show.re.var= TRUE, 
                  pred.labels =c("(Intercept)","Compact", "Marine Debris", "Rake", "Temperature"),
                  dv.labels= "Treatment and Temperature Effects on Burrow Diameter")

effects_treatment<- effects::effect(term= "mean.temp", mod= count.mod.final)
summary(effects_treatment)

df_treat <- as.data.frame(effects_treatment)

#Using predict()

temp_d <- expand.grid(
  mean.temp = seq(10,25,0.1),
  trt = d.comp3$trt, 
  site = d.comp3$site
)
temp_d

temp_d$count.pred <-predict(count.mod.final, newdata = temp_d)
temp_d$dia.pred <-predict(dia.model.final, newdata = temp_d)

count.pred <- ggplot(data=temp_d, aes(x = mean.temp, y = count.pred, colour = trt)) +
  geom_jitter(data=d.comp3, aes(x = mean.temp, y = count.diff, color=trt), width=0.05, alpha=0.5) +
  geom_line(size=1, alpha=0.5)+
  geom_hline(yintercept=0, size=1, alpha=0.5) +
  labs(x="Average Air Temperature (°C)", y= "Burrow Density Compared to First Day") +
  guides(color=guide_legend(title="Treatment"))+
  scale_colour_manual(labels = c("Control", "Compact", "Marine Debris", "Rake"), values = c("#999999", "#56B4E9","#009E73", "#D55E00")) +
  theme(text = element_text(size = 13))  

count.pred

dia.pred <- ggplot(data=temp_d, aes(x = mean.temp, y = dia.pred, colour = trt)) +
  geom_jitter(data=d.comp3, aes(x = mean.temp, y = dia.diff, color=trt), width=0.05, alpha=0.5) +
  geom_line(size=1, alpha=0.5)+
  geom_hline(yintercept=0, size=1, alpha=0.5) +
  labs(x="Average Air Temperature (°C)", y= "Burrow Diameter Compared to First Day") +
  guides(color=guide_legend(title="Treatment"))+
  scale_colour_manual(labels = c("Control", "Compact", "Marine Debris", "Rake"), values = c("#999999", "#56B4E9","#009E73", "#D55E00")) +
  theme(text = element_text(size = 13))  

dia.pred

#Goofin - using the day before as a comp 

d.comp_day <- d.agg %>%
  # Convert the date column to a Date object if it's not already
  mutate(date = as.Date(date)) %>%
  # Arrange the data by date, site, and treatment
  arrange(date, site, trt) %>%
  # Group by site and treatment
  group_by(site, trt) %>%
  # Calculate the difference in counts between each date and the first date
  mutate(count.diff = count - lag(count)) %>%
  mutate(dia.diff = mean.burrow.diameter-first(mean.burrow.diameter))

d.comp2_day <- d.comp_day %>%
  arrange(date) %>%
  filter(row_number() > 1)

write.csv(d.comp2_day, "d_comp2_day.csv", row.names=FALSE)
#I need to change the mean.temp NA's to the value from that day and I can't do that in R, so I wrote the .csv and did it manually 

d.comp3_day <- d_comp2_day %>%
  # Use mutate to replace a specific value in the column
  mutate(trt = ifelse(trt == "Control", "aControl", trt)) %>%
  mutate(trt = ifelse(trt == "Marine.Debris", "Marine_debris", trt))


cmod15 <- lm(count.diff~site, data=d.comp3_day)
summary(cmod1)
plot(cmod1, ask=F)
#Variance of burrow count comp appears homogeneous among replicates 

cmod16 <- lm(count.diff~site+mean.temp, data=d.comp3_day)
anova(cmod15, cmod16)

cmod17 <-lm(count.diff~trt + mean.temp + trt*mean.temp, data=d.comp3_day)
summary(cmod17)
plot(d.comp3_day$site, residuals(cmod17, type='pearson'))

cmod18 <- lme(count.diff~trt + mean.temp + trt*mean.temp, random= ~1|site, method = 'REML', data=d.comp3_day)
summary(cmod18)
plot(cmod18)

gls_cmod18 <- gls(count.diff~trt + mean.temp + trt*mean.temp, data=d.comp3_day)
anova(gls_cmod18, cmod18)

#No need for random effect with this analysis


#Refit model with ML and drop interaction 
cmod19 <- lm(count.diff~trt + mean.temp + trt*mean.temp, data=d.comp3_day)
cmod20 <- lm(count.diff~trt + mean.temp, data=d.comp3_day)
anova(cmod19,cmod20)
#So the interaction term can be dropped 

cmod19 <- lm(count.diff~trt + mean.temp,data=d.comp3_day)
cmod19_droptreat <- lm(count.diff~mean.temp, data=d.comp3_day)
cmod19_droptemp <- lm(count.diff~trt, data=d.comp3_day)

summary(cmod19)
summary(cmod19_droptreat)
summary(cmod19_droptemp)

anova(cmod19,cmod19_droptemp)
anova(cmod19,cmod19_droptreat)
#So treatment is not significant here when comparing to previous day 


#Just the last 3 days 
d.comp.goof <- d.comp3 %>%
  arrange(date) %>%
  filter(row_number() > 80)

count.goof <-lme(count.diff~trt, random= ~1|site, method = 'REML', data=d.comp.goof)


sjPlot::plot_model(count.goof, 
                   axis.labels=c("Rake", "Marine Debris", "Compact"),
                   show.values=TRUE, show.p=TRUE,
                   title="Effect of Treatments on Burrow Count")

temp_goof <- expand.grid(
  mean.temp = seq(10,25,0.1),
  trt = d.comp.goof$trt, 
  site = d.comp.goof$site
)
temp_goof

temp_goof$count.pred <-predict(count.goof, newdata = temp_goof)

ggplot(data=temp_goof, aes(x = mean.temp, y = count.pred, colour = trt)) +
  geom_jitter(data=d.comp3, aes(x = mean.temp, y = count.diff, color=trt), width=0.05, alpha=0.5) +
  geom_line(size=1, alpha=0.5)+
  geom_hline(yintercept=0, size=1, alpha=0.5) +
  labs(x="Average Air Temperature (°C)", y= "Burrow Density Compared to First Day") +
  guides(color=guide_legend(title="Treatment"))+
  scale_colour_manual(labels = c("Control", "Compact", "Marine Debris", "Rake"), values = c("#999999", "#56B4E9","#009E73", "#D55E00")) +
  theme(text = element_text(size = 13))  
