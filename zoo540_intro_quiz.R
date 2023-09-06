#Just plotting grouse at different wind levels
plot(data = grouse_data, GROUSE~WIND)

#Change up how many graphs are displayed 
par(mfrow = c(2, 2))

#Diagnostics 
plot(lm(data=grouse_data, GROUSE~WIND))

par(mfrow = c(1, 1))
#Make it an object
grouse_reg <- lm(data=grouse_data, GROUSE~WIND)
grouse_reg

plot(data=grouse_data, GROUSE~WIND)
abline(grouse_reg)

