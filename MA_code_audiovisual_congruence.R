#Main q: are children regardless of anything sensitive to congruent audio-visual 
#information? How big is the effect? Publication bias, etc.

#SD = squareroot(N) * (upperCI - lowerCI / 3.92), i.e 2 x 1.96 standard errors wide.
#SE = SD / squareroot(n) ; SD = SE * squareroot(n)

#convert t-statistic to effect size:
#for one-sample t:
  #Cohen's d = t/Sqrt(n)

library('tidyverse')
library('dplyr')

#import data
MA_data <- read.csv('/Users/au620441/Documents/GitHub/MyGitRepoCCox1/MA_audiovisual_congruence/MA_audiovisual_congruence/MA_Data_Audiovisual_Congruence4.csv')
MA_data <- as_tibble(MA_data)

#calculate average SD_1 to use for missing values:
MA_data %>%
  filter(SD_1 < 1) %>%
  summarise(mean(SD_1,na.rm=TRUE))

#use mean and sd to calculate d:
MA_data <- MA_data %>%
  mutate(MA_data, cohen_d = (x_1-x_2) / SD_1) %>%
  relocate(cohen_d, .after = t)

#use t-values to calculate d:
MA_data <- MA_data %>%
  mutate(MA_data, cohen_d2 = t / sqrt(n_1)) %>%
  relocate(cohen_d2, .after = cohen_d)

