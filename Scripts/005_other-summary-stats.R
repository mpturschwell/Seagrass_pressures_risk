library(dplyr)

dat2000 <- read.csv("Data/Data_seagrass_trends_pressures_2000-2010.csv")
#--------------------------------------------------------#
# Calculate frequency of each trajectory by bioregion 
#--------------------------------------------------------#
freq.table <-dat2000 %>% 
  group_by(bioregion) %>%
  count(category) %>%
  mutate(freq = n / sum(n))
freq.table

#--------------------------------------------------------#
# Count extreme values > 100% change per decade
#--------------------------------------------------------#
highdata <- dat2000[ which(dat2000$IRtrend_X50. > 0.1), ]
lowdata <- dat2000[ which(dat2000$IRtrend_X50. < -0.1), ]

