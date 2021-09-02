library(dplyr)

dat2000 <- read.csv("Data/Data_seagrass_trends_pressures_2000-2010.csv")

# Calculate frequency of each trajectory by bioregion 
freq.table <-dat2000 %>% 
  group_by(bioregion) %>%
  count(category) %>%
  mutate(freq = n / sum(n))
freq.table
