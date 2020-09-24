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
MA_data <- read.csv('/Users/au620441/Documents/GitHub/MyGitRepoCCox1/MA_audiovisual_congruence/MA_audiovisual_congruence/MA_Data_Audiovisual_Congruence6.csv')
MA_data <- as_tibble(MA_data)

#use mean and sd to calculate d:
MA_data <- MA_data %>%
  #calculate cohen's d with mean + sds:
  mutate(MA_data, cohen_d_mean_sd = (x_1-x_2) / SD_1) %>%
  #calculate cohen's d with t values:
  mutate(MA_data, cohen_d_t = t / sqrt(n_1)) %>%
  #merge the two above columns and include d, so there are no missing values:
  mutate(MA_data, cohen_d_final = coalesce(cohen_d_mean_sd,cohen_d_t, d)) %>%
  #position the columns for comparison:
  relocate(cohen_d_mean_sd, .before = study_ID) %>%
  relocate(cohen_d_t, .after = cohen_d_mean_sd) %>%
  relocate(cohen_d_final, .after = cohen_d_t)

#code for meta-analysis:
library(metafor)






