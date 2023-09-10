# This problem set introduces you to some of the basic mechanics of R. It is based on Peng (2016) R Programming for Data Science; Chapt headings below correspond to chapter in this book. You can download the book at  http://leanpub.com/rprogramming.

# An additional source that might be helpful is https://environmentalcomputing.net/about-this-site/. This presents similar material but in a different order. Also, it does more in dplyr, tidyr, and ggplot.

# PS1 is in two parts, PS1a and PS1b. The first is due by midnight on Sunday, 10 September. For this, I'd like all individuals to turn in their R own code. Part PS1b is due by midnight on Sunday, 17 September, and PS1b can be done in groups. Even though you do it in groups, each individual should submit R code on Canvas. Below I've pasted PS1a and PS1b together. For PS1a you can submit the entire file for PS1a (with the PS1b section left blank). For PS1b you can leave in your PS1a.

# For the homework, ONLY TURN IN THE R CODE THAT YOU USED. Start with the code from PS1_R_Programming_2023.R and add to it anything you need. If you are adding to my code, identify you new code (so I can find it easily) by placing it between marker rows #~~~~~~~~~~~~~~~~~~~~~~~~~~~. There is no need for you to write down detailed answers, since we will discuss everything in class. Do write down questions that you have, though.

# As a departure from the norm for later problem sets, EVERYBODY SHOULD TURN IN A COPY OF THEIR R CODE. I know that I will likely get replicate version from those of you i the same group. Still, I want to make sure everybody goes through all of the steps of submitting R code.

# The headings below refer to the corresponding chapter in Peng (2016); even though not all chapters are included, you should read all chapters.

# Finally, some of these questions might seem pointless. They are just to get you familiar with the tools you'll need in class. Also, some questions ask you to do the same thing in two or more ways. There are many ways to do the same simple things in R, and knowing different ways will help you work with other people in R.


# Due: 10Sep23 (PS1a) and 17Sep23 (PS1b)

#################################################################
#################################################################
# start PS1a
#################################################################
#################################################################


#################################################################
# Chapt 5. Nuts and Bolts
#################################################################

#1. Create a 2 x 5 matrix containing 1:10, with 1:5 on the first row and 6:10 on the second row. (Bonus: Come up with two or more ways to do this.)
matrix_1 <- matrix(1:10, nrow=2, ncol=5)
matrix_1
#This reads vertically, which is not what we want! Instead, let's try this:
x1 <- 1:5
y1 <- 6:10
rbind(x1,y1)
#That's better!

#2. Label the rows of your matrix "A" and "B", and the columns "a", ..."e".
#Let's make an object 
matrix_2 <- rbind(x1,y1)
matrix_2
#Now let's change the row and column names 
rownames(matrix_2) <- c("A","B")
colnames(matrix_2) <- letters[1:5]
#Checking it out 
matrix_2

#3. Convert the values of your matrix to double.
matrix_2*2

#4. Create a data.frame named "df" with a column named "site" containing values 1:5 as factors and a column named "value" containing values 100, 200, ..., 500. Then give it row.names "a", ..., "e".
df <- data.frame(site = 1:5, value = c(100,200,300,400,500))
df
rownames(df) <- letters[1:5]
df

#Question: is there a way to more easily assign names a:e where I don't have to type out every row name?
#There is indeed! We can use the argument 'letters[x:y]' with how many letters you want! 

#5. Rename the variables (columns) in df "site.cat" and "x".
colnames(df) <- c("site.cat", "x")
df

#6. Convert df to a list and extract the first element of the list by number and then by name ("site.cat").
list_1 <- as.list(df)
list_1

first_element <- list_1$site.cat
second_element <- list_1$x

first_element
second_element
#Is this what was meant by extract? Create new objects for each variable & values?

#################################################################
# Chapt 6. Uploading data
#################################################################

# Metadata for "grouse_data.csv"
# These data are simulated to have similar characteristics as the original, real data that will be discussed in class for PS2. I thought I'd use it as an example here, although you don't need to know any details about it.

# ROUTE
# IDs for 50 roadside routes. 

# STATION
# IDs for each survey station, with up to 8 STATIONS per ROUTE.

# LAT
# X coordinate of survey station. UTM zone 15N NAD83 datum.

# LONG
# Y coordinate of survey station. UTM zone 15N NAD83 datum.

# WIND
# Wind speed (km/hour) recorded at 1.4m above ground at the end of each observation.

# TEMP
# Temperature (°C) recorded at 1.4m above ground at the end of each observation.

# GROUSE
# Detection/non-detection of Ruffed Grouse (1 = detected, 0 = not detected).

#7. Create a data.frame called d by reading data from file = "grouse_data.csv". Make the variables "ROUTE" and "STATION" factors, and "GROUSE" integers
#Tony - 
#I use github to host my different R projects, so reading in data may be a little different! 
#I like this better as a way to organize my working directories, but can change it for this class if needed! 
library(readr)
d <- read_csv("grouse_data.csv")
#view(d)
class(d$ROUTE)
#let's fix this to be a factor 
d$ROUTE <- as.factor(d$ROUTE)
class(d$ROUTE)
#Much better, onto the next 
class(d$STATION)
#Same change here 
d$STATION <- as.factor(d$STATION)
class(d$STATION)
#Finally, change grouse 
class(d$GROUSE)
d$GROUSE <- as.integer(d$GROUSE)
class(d$GROUSE)
#################################################################
# Chapt 10. Subsetting data
#################################################################

#8. Select the subset of d with (a) ROUTE == 1, (b) LAT greater than mean(LAT), (c) LAT greater than median(LAT). Now compute the number of locations (rows) in d with LAT greater than 570000.
#For (a)
subset(d, ROUTE == 1)
#For (b)
subset(d, LAT > mean(LAT))
#For (c)
subset(d, LAT > median(LAT))
#Counting rows with LAT > 570000
sum(d$LAT > 570000)
#There are 291 LAT values > 570000
#I thought there was a count function as well?

#9. Create a new variable "NORTHERN" that is TRUE for the more northerly half of the stations and FALSE otherwise.
#How do we define northern half and southern half? 
#Do we just use the median LAT to designate N vs. S? 
median(d$LAT)
#The median LAT is 578273, so this will be our cutoff 
d$NORTHERN <- if(d$LAT > median(d$LAT)){
  print("TRUE")
} else {
    print("FALSE")
}
#Well R doesn't like this very much
library(tidyverse)
#Reverting to tidyverse - sorry! 
d <- d %>%
  add_column(NORTHERN = 
               if_else(.$LAT > median(.$LAT), TRUE, FALSE),
             .after = "GROUSE")

#10. How many different routes are there in d?
length(unique(d$ROUTE))
#There's 50 different routes in d! 

#11. Construct a matrix d.matrix that has the same information as d. Compute the number of routes using d.matrix and the unique() function.
d.matrix <- read_csv("grouse_data.csv")
d.matrix
length(unique(d.matrix$ROUTE))
#There are 50 unique routes! 

#################################################################
# Chapt 11. Vectorized Operations
#################################################################

#12. Create a new variable "QUADRANT" that divides the stations into four groups, with values 0 for the northwest group, 1 for the northeast group, 2 for the southwest group, and 3 for the southeast group.
#Wait, does Longitude increase east to west or west to east? 
#Maybe I need an intro to geography class 
#West is more, east is less (I think in the US)
d <- d %>%
  add_column(QUADRANT = 
               case_when(.$LAT > median(.$LAT) & .$LONG > median(.$LONG) ~ "0", 
                         .$LAT > median(.$LAT) & .$LONG < median(.$LONG) ~ "1",
                         .$LAT < median(.$LAT) & .$LONG > median(.$LONG) ~ "2",
                         .$LAT < median(.$LAT) & .$LONG < median(.$LONG) ~ "3"),
             .after = "NORTHERN")
#Watch your commas here! Took me too long to realize i need a ~ and not ,

#################################################################
# Chapt 12. Dates and times
#################################################################

# You should read this chapter to be aware these features exist.

#################################################################
# Chapt 13. Managing data.frames
#################################################################
#dplyr should already be loaded after libbing tidyverse! 

#13. Select columns of d that start with "L" (i.e., LAT and LONG). Do this both using and not using dplyr.
#No dplyr 
d[,startsWith(colnames(d), "L")]
#Dplyr 
d %>%
  select(starts_with("L"))

#14. Select rows of d that contain the highest 50% of the values of LAT. Do this both using and not using dplyr. Do you see any differences in the data.frames that are produced (check the row.names)?
subset(d, d$LAT > mean(d$LAT))

dplyr::filter(d, d$LAT > mean(d$LAT))
#No differences that I can see? 

#15. Select rows of d that contain ROUTE > 45. This is much trickier than you might think, because ROUTE is a factor!
#Alrighty, do this without turning route back into numeric class 

dplyr::filter(d, ROUTE %in% c("45", "46", "47", "48", "49", "50"))
#Long but it works! 

#16. Reorder rows by LAT to create a new data.frame d.LAT. Then reorder d.LAT by ROUTE
d.LAT <- d %>%
  arrange(LAT)
d.LAT
d.LAT <- d.LAT %>%
  arrange(LAT, ROUTE)
d.LAT
#I think this is what was meant? Increasing LAT, but within same LATS have increasing ROUTE?

#17. Rename the column in d from ROUTE.num to num.ROUTE. Do this both using and not using dplyr.
#Did we ever rename this above? 
#No dplyr
colnames(d)[1] ="num.ROUTE"
#Dplyr 
d <- d %>%
  rename("ROUTE.num" = "num.ROUTE")
#Important to note that the new column name goes first

# NOTE: I am assuming that you carried through all the variables in data.frame d from all the questions above. If not, just rerun all your code up to this point.

#18. Add a variable to d called exp.GROUSE than contains exp(GROUSE). Do this both using and not using dplyr.
#NO dplyr
d$exp.GROUSE <- exp(d$GROUSE)
#Dplyr 
d <- d %>%
  add_column(dexp.GROUSE = exp(d$GROUSE))
#Renamed this to dexp.GROUSE to make sure it was correctly calculating and adding the column 

#19. Create a data.frame called d.ROUTE that contains the mean value of GROUSE for each ROUTE. Do this both using and not using dplyr. Bonus: create a new data.frame d.NA which is the same as d but with STATION 1 in ROUTE 1 having a value of NA (i.e., d$ROUTE[1] <- NA). Then perform the same operation as you just did for d. Note that you hve to figure out a way to get mean(GROUSE) for ROUTE == 1 factoring out the NA at STATION == 1.
#Wait, isn't grouse binary? 
#No dplyr 
d.ROUTE <- aggregate(d$GROUSE, by=list(d$ROUTE.num), FUN=mean)
#This is a little gross 
#Using dplyr 
dd.ROUTE <- d %>%
  group_by(ROUTE.num) %>%
  summarize(GROUSE=mean(GROUSE), 
  .groups='drop')

#20. Perform the same operation as in #19 using piping in dplyr.
#Is this what I did above? How should I do this not using piping? 

#################################################################
#################################################################
# end PS1a, start PS1b
#################################################################
#################################################################

#################################################################
# Chapt 14. Control Structures
#################################################################
# The tasks requested below could/should be done as in the previous questions. Here, I'm using the same tasks to illustrate the basics of control structures. The basic structures I'm illustrating are very common, and you will see them repeatedly throughout the course. This also shows that you can pretty much do anything with loops, which isn't surprising given that the computer code underlying functions like group_by() and aggregate() uses loops.

#21. Create a data.frame called d.ROUTE that contains the mean value of GROUSE for each ROUTE using a for loop (i.e., group_by() and aggregate() aren't allowed). This involves first creating a data.frame and then looping through the values of ROUTE to fill in values. I've given you the first line creating the data.frame
d.ROUTE <- data.frame(ROUTE = unique(d$num.ROUTE))

#22. Create a data.frame called d.ROUTE that contains the mean value of GROUSE for each ROUTE using a while loop (i.e., group_by() and aggregate() aren't allowed). 
d.ROUTE <- data.frame(ROUTE = unique(d$num.ROUTE))

#23. Create a data.frame called d.ROUTE.NA from d.NA that contains the mean value of GROUSE for each ROUTE using a for loop, in which the ROUTE with an NA is given the value NaN. You will need the is.na() function for this.
d.ROUTE.NA <- data.frame(ROUTE = unique(d.NA$num.ROUTE))

#################################################################
# Chapt 15. Functions
#################################################################

#24. Write a function locate_station() that returns the LAT and LONG for a ROUTE and STATION in d. Note that you are best to use num.ROUTE, not ROUTE.

#25. Write a function distance_station() that returns a new data.frame containing the Euclidean distance between a specified ROUTE and STATION in d and all other routes and stations. Don't bother about converting to meters -- just calculate distance in terms of UTM values (LAT and LONG). Note that you can use the function locate_station() in your new function.

# 26. Write a function plot_station() that plots the location of stations, with the size of the point for each station decreasing with increasing distance from a focal station specified as input to the function.

#################################################################
# Chapt 16. Scoping rules in R
#################################################################

# You should read this chapter to be warned. However, I'm going to try to make sure all the code you use in class is properly scoped. The bottom line is that it is often best to define variables in the definition of a function.

# 27. Take a look at the following two sets of code. BEFORE running the code, predict what the outcomes of f1 and f2 will be.

y <- 10
f1 <- function(x) {
	y <- 2
	y^2 + g1(x)
}
g1 <- function(x) {
	x*y
}
f1(3)

y <- 10
f2 <- function(x) {
	g2 <- function(x) {
		x*y
	}
	y <- 2
	y^2 + g2(x)
}
f2(3)

#################################################################
# Chapt 17. Standards
#################################################################

# These are good principles (but don't expect me to always follow them).

#################################################################
# Chapt 18. Loop functions
#################################################################

# 28. Use tapply() to create the same date.frame as aggregate(GROUSE ~ num.ROUTE, data=d, FUN=mean). You can also try this using ROUTE rather than num.ROUTE, and you will see some problems.
aggregate(GROUSE ~ num.ROUTE, data=d, FUN=mean)

# 29. First, add columns to data.frame d which give the distance of each station from two locations (e.g., ROUTE = 1, STATION = 1; and ROUTE = 40, STATION = 1). Use apply() to create an additional column that is the sum of these distances.

# 30. Use lapply() and split() to create the same date.frame as aggregate(GROUSE ~ num.ROUTE, data=d, FUN=mean). You can also try this using ROUTE rather than num.ROUTE, and you will see some problems.
aggregate(GROUSE ~ num.ROUTE, data=d, FUN=mean)

#################################################################
# Chapts 19/20/21. Regular expressions, Debugging and Profiling
#################################################################

# You should read these chapters to be aware these features exist.

#################################################################
# Chapt 22. Simulation
#################################################################

# 31. Plot a histogram of 1000 values simulated from a Gaussian distribution with mean = 0 and standard deviation = 1. Then plot the probability density function to compare them. You will need the hist() function (with the freq=F option).

# 32. Plot a histogram of 1000 values simulated from a Poisson distribution with lambda = 3. Then plot the probability density (or mass) function to compare them. You will need the hist() function (with the freq=F option).

