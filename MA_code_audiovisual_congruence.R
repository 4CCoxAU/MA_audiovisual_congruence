library('tidyverse')
library('dplyr')
library('moments')
library('car')
library('metafor')
library('esc')
library('ggplot2')
library('ggridges')
library('brms')

#import data:
MA_data <- read.csv('/Users/au620441/Documents/GitHub/MyGitRepoCCox1/MA_audiovisual_congruence/MA_audiovisual_congruence/MA_Data_Audiovisual_Congruence_final_subset_2.csv')
MA_data <- as_tibble(MA_data)

#specifications of priors for model:
priors <- c(prior(normal(0,1), class = Intercept),
            prior(cauchy(0,0.5), class = sd))

#model with fixed and random effects:
m.brm <- brm(hedge_g|se(se_hedge_g) ~ 1 + mean_age_1 + test_lang + (1|study_ID),
             data = MA_data,
             prior = priors,
             iter = 10000, 
             warmup = 2000)

#checks & outcomes:
pp_check(m.brm)
mcmc_plot(m.brm, type='hist')
summary(m.brm)
ranef(m.brm)

#plot posterior samples for effect size:
post.samples <- posterior_samples(m.brm, c("^b", "^sd"))
names(post.samples) <- c("smd", "tau")

ggplot(aes(x = smd), data = post.samples) +
  geom_density(fill = "lightblue", color = "lightblue", alpha = 0.7) +
  geom_point(y = 0, x = mean(post.samples$smd)) +
  labs(x = expression(italic(SMD)),
       y = element_blank()) +
  theme_minimal()

#probability that effect size is under 0.3:
smd.ecdf <- ecdf(post.samples$smd)
smd.ecdf(0.3)


#use esc to create forest plot for random effects:
full.model <- rma.mv(hedge_g, 
                     se_hedge_g,
                     random = ~ 1 | study_ID, 
                     tdist = TRUE, 
                     data = MA_data,
                     method = "REML")
summary(full.model)
forest(full.model)
funnel(full.model)


