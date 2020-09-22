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
MA_data <- read.csv('/Users/au620441/Documents/GitHub/MyGitRepoCCox1/MA_audiovisual_congruence/MA_audiovisual_congruence/MA_Data_Audiovisual_Congruence3.csv')
MA_data <- as_tibble(MA_data)

#use t-values to calculate d:
MA_data <- MA_data %>%
  mutate(MA_data, cohen_d = t / sqrt(n_1)) %>%
  relocate(cohen_d, .after = t)

#calculate average SD_1 to use for missing values:
MA_data %>%
  filter(SD_1 < 1) %>%
  summarise(mean(SD_1,na.rm=TRUE))

#calculate average SD_2 (baseline condition) to use for missing values:
MA_data %>%
  filter(SD_2 < 1) %>%
  summarise(mean(SD_2,na.rm=TRUE))
    


