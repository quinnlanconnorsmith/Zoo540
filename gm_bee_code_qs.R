library(tidyverse)
data <- read.csv('bees_combined.csv') %>% filter(!grepl("date", date)) %>% mutate(trt =  ifelse(grepl("2|3", cameraID), "ash", "control"))


#Make an integrated timestamp
start.date <- parse_date_time('2023-09-06 06:24:28', "%Y-%m-%d %H:%M:%S")
data$timestamp <- paste(data$date, substr(data$short, 13, 18), sep = "_") 
data$timestamp <- parse_date_time(data$timestamp, "ymdHMS")
data$time.num <- as.numeric(difftime(data$timestamp, start.date, units = 'days'))

#conf and prob1 are the detection thresholds, the experiment didn’t start until the 9th even though data collection starts on the 6th: 
  det_thresh <- 0.7
class_thresh <- 0.7

#high quality data (detection and class conf > 0.7)
hq.data <- data %>% filter(conf > det_thresh & prob1 > class_thresh, 
                           timestamp > "2023-09-09 11:00:00")
#Q code 

# Convert the timestamp column to a POSIXct object
hq.data$timestamp <- as.POSIXct(hq.data$timestamp, format="%Y-%m-%d %H:%M:%S")

# Create a new column for 15-minute bins
qdata <- hq.data %>%
  mutate(bin = cut(timestamp, breaks = "15 min"))

# Group by date, cameraID, and the 15-minute bin
qdata2 <- qdata %>%
  group_by(date = date(timestamp), cameraID, bin) %>%
  summarise(count = n()) %>%
  ungroup()
