library(tidyverse)
library(dplyr)
library(moments)
library(metafor)
library(ggplot2)
library(ggridges)
library(brms)
library(metaviz)
library(mice)
library(lme4)
library(pwr)
library(brmstools)
library(PublicationBias)
library(loo)
library(rstan)

#import data:
MA_data <- read.csv('/Users/au620441/Documents/GitHub/MyGitRepoCCox1/MA_audiovisual_congruence/MA_audiovisual_congruence/Final_Data/Redo_Analysis_12-10/hedge_es_final_calculation.csv')
MA_data <- as_tibble(MA_data)

#Calculations for missing data:

#Overview of missing data:
md.pattern(MA_data)

#first mice is just to get the parameters, so these can be restricted to positive values:
MA_data_imp <- mice(MA_data,m=1,maxit=25,meth='norm',seed=7, print=F) #norm is Bayesian linear regression.
post <- MA_data_imp$post
meth <- MA_data_imp$meth
post["se_hedge_g"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0.1445317, 1))"

#post-process them, so they are restricted to positive range
MA_data_imp <- mice(MA_data, meth=meth, post=post, print=FALSE, m=20, maxit=25)

#check values to make sure they all make sense
MA_data_imp$imp$se_hedge_g
densityplot(MA_data_imp)

#Bayesian hierarchical model of effect size, using imputed dataset:

#Model 1:
priors1 <- c(prior(normal(0, 0.5), class = Intercept),
             prior(gamma(2, 0.1), class = nu))
brm.student_baseline <- 
  brm_multiple(data = MA_data_imp, family = student,
               hedge_g|se(se_hedge_g) ~ 1 + (1|study_ID/expt_condition),
               prior = priors1,
               iter = 20000, 
               warmup = 2000,
               cores=4,
               control = list(adapt_delta = 0.99))
pp_check(brm.student_baseline)
summary(brm.student_baseline)

#Model 2:
priors2 <- c(prior(normal(0, 0.5), class = Intercept),
             prior(normal(0, 0.5), class = b),
             prior(gamma(2, 0.1), class = nu))
brm.student_age <- 
  brm_multiple(data = MA_data_imp, family = student,
               hedge_g|se(se_hedge_g) ~ 1 + mean_age_1 + (1|study_ID/expt_condition),
               prior = priors2,
               iter = 20000, 
               warmup = 2000,
               cores=4,
               control = list(adapt_delta = 0.99))
pp_check(brm.student_age)
plot(conditional_effects(brm.student_age), points = TRUE)
summary(brm.student_age)

#Model 3:
brm.student_lang <- 
  brm_multiple(data = MA_data_imp, family = student,
               hedge_g|se(se_hedge_g) ~ 1 + test_lang + (1|study_ID/expt_condition),
               prior = priors2,
               iter = 20000, 
               warmup = 2000,
               cores=4,
               control = list(adapt_delta = 0.99))
pp_check(brm.student_lang)
plot(conditional_effects(brm.student_lang), points = TRUE)
summary(brm.student_lang)

#Model 4:
brm.student_stimuli <- 
  brm_multiple(data = MA_data_imp, family = student,
               hedge_g|se(se_hedge_g) ~ 1 + stimuli_type + (1|study_ID/expt_condition),
               prior = priors2,
               iter = 20000, 
               warmup = 2000,
               cores=4,
               control = list(adapt_delta = 0.99))
pp_check(brm.student_stimuli)
plot(conditional_effects(brm.student_stimuli), points = TRUE)
summary(brm.student_stimuli)

#Model 5:
brm.student_interaction <- 
  brm_multiple(data = MA_data_imp, family = student,
               hedge_g|se(se_hedge_g) ~ 1 + mean_age_1*test_lang + (1|study_ID/expt_condition),
               prior = priors2,
               iter = 20000, 
               warmup = 2000,
               cores=4,
               control = list(adapt_delta = 0.99))
pp_check(brm.student_interaction)
plot(conditional_effects(brm.student_interaction), points = TRUE)
summary(brm.student_interaction)

#plot the influence of the moderators:
plot(conditional_effects(brm.student_age, spaghetti = T, nsamples = 250), points = T, point_args = c(alpha = 0.9, size = 1.7), mean=F)
plot(conditional_effects(brm.student_lang), title = 'DSDSDSD')
plot(conditional_effects(brm.student_stimuli))
plot(conditional_effects(brm.student_interaction,
                         effects = "mean_age_1:test_lang", 
                         spaghetti = T, nsamples = 300),points = T, point_args = c(alpha = 0.9, size = 2), mean = F)

#Check fits and model comparison:
brm.student_baseline <- add_ic(brm.student_baseline, ic="loo")
brm.student_age <- add_ic(brm.student_age, ic="loo")
brm.student_lang <- add_ic(brm.student_lang, ic="loo")
brm.student_stimuli <- add_ic(brm.student_stimuli, ic="loo")
brm.student_interaction <- add_ic(brm.student_interaction, ic="loo")

loo_model2 <- loo(brm.student_lang, moment_match = T)
plot(loo_model2, label_points = T)
pareto_k_influence_values(loo_model)

loo_model_weights(brm.student_baseline, brm.student_age)
loo_model_weights(brm.student_baseline, brm.student_lang)
loo_model_weights(brm.student_baseline, brm.student_stimuli)
loo_model_weights(brm.student_baseline, brm.student_interaction)

#save models:
saveRDS(brm.student_baseline, "brm.student_baseline.rds")
saveRDS(brm.student_age, "brm.student_age.rds")
saveRDS(brm.student_lang, "brm.student_lang.rds")
saveRDS(brm.student_stimuli, "brm.student_stimuli.rds")
saveRDS(brm.student_interaction, "brm.student_interaction.rds")

#Create forest plot:
forest <- forest(
  model = brm.student_baseline, 
  sort = TRUE,
  fill_ridge = "lightsteelblue4")
forest + 
  scale_x_continuous("Hedges' g", limits = c(-0.75, 2.75)) +
  ylab(' ') +
  ggtitle("Forest Plot of Estimated Effect Sizes") +
  theme_forest()

#Examine publication bias: 

#create average of 100 imputations: use only for significance funnel plot & publication bias analysis
MA_data_imp_for_average <- mice(MA_data, meth=meth, post=post, print=FALSE, m=100, maxit=25)
all_data_imputations <- as.tibble(MA_data_imp_for_average$imp$se_hedge_g)
#mean across 100 data imputations:
mean_mi <- rowMeans(average_mi)
mean_mi <- as.tibble(mean_mi)
#combine original with imputed sd values:
original_values <- na.omit(MA_data$se_hedge_g)
MA_data_average_imp <- MA_data %>%
  mutate(se_hedge_g = c(original_values, mean_mi$value)) %>%
  relocate(se_hedge_g, .before = n_1)

#make sensitivity plot, as in Mathur & VanderWeele (2020):
eta.list = as.list( c( 200, 150, 100, 50, 40, 30, 20, rev( seq(1,15,1) ) ) )
res.list = lapply( eta.list, function(x) {
  cat("\n Working on eta = ", x)
  return( corrected_meta( yi = MA_data_average_imp$hedge_g,
                          vi = MA_data_average_imp$se_hedge_g,
                          eta = x,
                          model = "robust",
                          favor.positive = TRUE,
                          clustervar = MA_data_average_imp$study_ID) )
}
)
# put results for each eta in a dataframe and plot:
res.df = as.data.frame( do.call( "rbind", res.list ) )
#plot the results:
ggplot( data = res.df, aes( x = eta, y = est ) ) + 
  geom_ribbon( data = res.df, aes( x = eta, ymin = lo, ymax = hi ), fill = "gray" ) +
  geom_line( lwd = 1.2 ) +
  xlab( bquote( eta ) ) +
  ylab( bquote( hat(mu)[eta] ) ) +
  theme_classic()


#create a significance funnel plot to examine publication bias:
sig_fun <- significance_funnel( yi = MA_data_average_imp$hedge_g,
                     vi = MA_data_average_imp$se_hedge_g,
                     xmin = min(MA_data_average_imp$hedge_g),
                     xmax = max(MA_data_average_imp$hedge_g),
                     ymin = 0.2,
                     ymax = max(sqrt(MA_data_average_imp$se_hedge_g)),
                     xlab = "Point Estimate of Effect Size",
                     ylab = "Standard Error of Effect Size",
                     favor.positive = T,
                     alpha.select = 0.05)
sig_fun +
  ggtitle("Significance Funnel Plot of Meta-Analytic Studies") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=14)) +
  geom_vline(xintercept=0.36, linetype="solid") +
  geom_vline(xintercept=0.238, linetype="dashed")


#calculate severity of publication bias needed to "explain away" the results:
svalue( yi = MA_data_average_imp$hedge_g,
        vi = MA_data_average_imp$se_hedge_g,
        q=0,
        clustervar = MA_data_average_imp$study_ID,
        model = "robust",
        alpha.select = 0.05,
        eta.grid.hi = 5,
        favor.positive = TRUE,
        CI.level = 0.95,
        small = TRUE,
        return.worst.meta = TRUE)


#power calculations:
pwr.t.test(n = , d = 0.33, sig.level = 0.05, power = 0.80, type = "one.sample")

#plot of posterior samples for effect size:
post.samples <- posterior_samples(m.brm, c("^b", "^sd"))
names(post.samples) <- c("smd", "tau")
mean <- mean(post.samples$smd)
summary(m.brm)

#plot the data:
ggplot(aes(x = smd), data = post.samples) +
  geom_density(fill = "lightsteelblue4", color = "black", alpha = 0.9) +
  geom_point(y = 0, x = mean(post.samples$smd)) +
  geom_vline(xintercept=mean) +
  geom_vline(xintercept=0.18, linetype="dotted") +
  geom_vline(xintercept=0.51, linetype="dotted") +
  ggtitle("Posterior Distribution for Pooled Effect Size Estimate") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=14)) +
  xlab('Effect Size') +
  ylab('Density')

ggplot(aes(x = tau), data = post.samples) +
  geom_density(fill = "lightblue", color = "lightblue", alpha = 0.7) +
  geom_point(y = 0, x = mean(post.samples$tau)) +
  labs(x = expression(italic(Tau)),
       y = element_blank()) +
  theme_minimal()

#probability that effect size is under 0.2:
smd.ecdf <- ecdf(post.samples$smd)
smd.ecdf(0.2)

#plot of age of infant participants
p <- MA_data_2 %>%
  ggplot( aes(x=mean_age_1)) +
  geom_histogram( binwidth=100, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("Age of Infant Participants in Meta-analytic Studies") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=14)) +
  xlab('Mean Age in Days') +
  ylab('Count')






#Below are deleted scenes:

"
#model to estimate se of effect size in studies:
m.brm_sd <- brm(se_hedge_g ~ n_1*hedge_g,
             data = MA_data,
             prior = priors,
             iter = 10000, 
             warmup = 2000,
             cores=4,
             control = list(adapt_delta = 0.80))

Posterior <- posterior_samples(m.brm_sd)
Posterior
Posterior$b_Intercept

data$Predicted_SE1 <- rnorm(1, mean(Posterior$b_Intercept), sd(Posterior$b_Intercept)) + rnorm(1, mean(Posterior$b_n_1, sd(Posterior$b_n_1))) + Posterior$b_hedge_g + Posterior$sigma
data$Predicted_SE2 <- rnorm(1, mean(Posterior$Intercept), sd(Posterior$Intercept)) + rnorm(1, mean(Posterior$b_n_1, sd(Posterior$b_n_1))  n_1 + b2  Hedges_g + b3  n_1  Hedges_g + sigma
predict(m, newdata = d_NA, allow_new_levels = TRUE, summary = TRUE))

#predict standard error based on model:
predicted_se <- predict(m.brm_sd)
predicted_se
dat <- as.data.frame(cbind(Y = standata(m.brm_sd)$Y, predicted_se))
ggplot(dat) + geom_point(aes(x = Estimate, y = Y))

#forest plot:
y <- data.frame(es=MA_data_2$hedge_g, se = MA_data_2$se_hedge_g, study_ID = MA_data_2$study_ID)
y <- y[order(MA_data_2$hedge_g),]
viz_forest(x = y[1:92, c("es", "se")], study_labels = y[1:92, c("study_ID")],summary_label = "Summary effect", xlab = "Hedges_g", variant = "rain", method = "DL")
"

#deleted scenes #2:
#frequentist model of data:
full.model <- rma.mv(hedge_g, 
                     se_hedge_g,
                     random = ~ 1 | study_ID/expt_condition, 
                     tdist = TRUE, 
                     data = MA_data_imp,
                     method = "REML")"
