pacman::p_load(
  #general
  tidyverse,
  here,
  dplyr,
  moments,
  glue,
  tidybayes,
  #metaanalysis
  metafor,
  metaviz,
  PublicationBias,
  #visualization
  ggplot2,
  ggridges,
  #modelfitting
  brms,
  lme4,
  mice,
  pwr)

#import data:
MA_data <- read.csv(here("hedge_es_final_calculation.csv"))
MA_data <- as_tibble(MA_data)

#Calculations for missing data:

#Overview of missing data:
md.pattern(MA_data)

#first mice is just to get the parameters, so these can be restricted to positive values:
MA_data_imp <- mice(MA_data,
                    m = 1,
                    maxit = 25,
                    meth = 'norm',
                    seed = 7, 
                    print = F) #norm is Bayesian linear regression.
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
baseline_f <- bf(hedge_g | se(se_hedge_g) ~ 1 + (1 | study_ID/expt_condition))

get_prior(baseline_f,
          data = MA_data, 
          family = student)

priors1 <- c(prior(normal(0, 0.5), class = Intercept),
             prior(normal(0, 0.2), class = sd),
             prior(gamma(2, 0.1), class = nu))

brm.student_baseline_prior <- 
  brm(
    baseline_f,
    data = subset(MA_data, !is.na(hedge_g)),
    family = student,
    prior = priors1,
    sample_prior = "only",
    file = "brm.student_baseline_prior",
    iter = 20000,
    warmup = 2000,
    cores = 2,
    chains = 2,
    control = list(adapt_delta = 0.99,
                   max_treedepth = 20))
pp_check(brm.student_baseline_prior, nsamples = 100)

brm.student_baseline <- 
  brm_multiple(
    hedge_g | se(se_hedge_g) ~ 1 + (1 | study_ID/expt_condition),
    save_pars = save_pars(all = TRUE),
    data = MA_data_imp, 
    family = student,
    prior = priors1,
    file = "brm.student_baseline",
    sample_prior = T,
    iter = 20000, 
    warmup = 2000,
    cores = 2,
    chains = 2,
    #backend = "cmdstanr",
    #threads = threading(2),
    control = list(
      adapt_delta = 0.99,
      max_treedepth = 20 ))
summary(brm.student_baseline)
pp_check(brm.student_baseline)

# Check posterior update (were the priors ok?)
Posterior <- posterior_samples(brm.student_baseline, pars = c(
  "prior_Intercept",
  "b_Intercept",
  "prior_sd_study_ID",
  "sd_study_ID__Intercept",
  "prior_sd_study_ID:expt_condition",
  "sd_study_ID:expt_condition__Intercept",
  "sigma",
  "prior_nu",
  "nu"
))
Posterior
ggplot(Posterior) +
  geom_density(aes(prior_Intercept), fill="blue", color="blue",alpha=0.3) +
  geom_density(aes(b_Intercept), fill="red", color="red",alpha=0.3) + 
  theme_classic()

ggplot(Posterior) +
  geom_density(aes(prior_nu), fill="blue", color="blue",alpha=0.3) +
  geom_density(aes(nu), fill="red", color="red",alpha=0.3) + 
  theme_classic()

ggplot(Posterior) +
  geom_density(aes(`prior_sd_study_ID:expt_condition`), fill="blue", color="blue",alpha=0.3) +
  geom_density(aes(`sd_study_ID:expt_condition__Intercept`), fill="red", color="red",alpha=0.3) + 
  theme_classic()

ggplot(Posterior) +
  geom_density(aes(`prior_sd_study_ID`), fill="blue", color="blue",alpha=0.3) +
  geom_density(aes(`sd_study_ID__Intercept`), fill="red", color="red",alpha=0.3) + 
  theme_classic()


summary(brm.student_baseline)

#Model 2:
priors2 <- c(prior(normal(0, 0.5), class = Intercept),
             prior(normal(0, 0.2), class = b),
             prior(normal(0, 0.2), class = sd),
             prior(gamma(2, 0.1), class = nu))

brm.student_age <- 
  brm_multiple(data = MA_data_imp, family = student,
               hedge_g|se(se_hedge_g) ~ 1 + mean_age_1 + (1|study_ID/expt_condition),
               prior = priors2,
               sample_prior = T,
               file = "brm.student_age",
               iter = 20000, 
               warmup = 2000,
               chains = 2,
               cores = 2,
               control = list(adapt_delta = 0.99))
pp_check(brm.student_age)

Posterior <- posterior_samples(brm.student_age, pars = c(
  "prior_Intercept",
  "b_Intercept",
  "prior_sd_study_ID",
  "sd_study_ID__Intercept",
  "prior_sd_study_ID:expt_condition",
  "sd_study_ID:expt_condition__Intercept",
  "sigma",
  "prior_nu",
  "nu"
))

pp_check(brm.student_age)

#model 2.5:
brm.student_age_mo <-
  brm_multiple(data = MA_data_imp, family = student,
               
               hedge_g|se(se_hedge_g) ~ 1 + mo(as.ordered(mean_age_1)) + (1|study_ID/expt_condition),
               prior = priors2,
               sample_prior = T,
               file = "brm.student_age_mo",
               iter = 20000,
               warmup = 2000,
               chains = 2,
               cores = 2,
               control = list(adapt_delta = 0.99))

pp_check(brm.student_age_mo, nsamples = 100)

brm.student_age <- add_criterion(brm.student_age, criterion = "loo")
brm.student_age_mo <- add_criterion(brm.student_age_mo, criterion = "loo")
loo_compare(brm.student_age,brm.student_age_mo)

Posterior <- posterior_samples(brm.student_age_mo, pars = c(
  "prior_Intercept",
  "b_Intercept",
  "prior_sd_study_ID",
  "sd_study_ID__Intercept",
  "prior_sd_study_ID:expt_condition",
  "sd_study_ID:expt_condition__Intercept",
  "sigma",
  "prior_nu",
  "nu"
))
#Model 3:
brm.student_lang <- 
  brm_multiple(data = MA_data_imp, family = student,
               hedge_g|se(se_hedge_g) ~ 1 + test_lang + (1|study_ID/expt_condition),
               prior = priors2,
               sample_prior = T,
               file = "brm.student_lang",
               iter = 20000, 
               warmup = 2000,
               chains = 2,
               cores = 2,
               control = list(adapt_delta = 0.99))
pp_check(brm.student_lang)
plot(conditional_effects(brm.student_lang), points = TRUE)
summary(brm.student_lang)

Posterior <- posterior_samples(brm.student_lang, pars = c(
  "prior_Intercept",
  "b_Intercept",
  "prior_sd_study_ID",
  "sd_study_ID__Intercept",
  "prior_sd_study_ID:expt_condition",
  "sd_study_ID:expt_condition__Intercept",
  "sigma",
  "prior_nu",
  "nu"
))


#Model 4:
brm.student_stimuli <- 
  brm_multiple(data = MA_data_imp, family = student,
               hedge_g|se(se_hedge_g) ~ 1 + stimuli_type + (1|study_ID/expt_condition),
               prior = priors2,
               sample_prior = T,
               file = "brm.student_stimuli",
               iter = 20000, 
               warmup = 2000,
               chains = 2,
               cores = 2,
               control = list(adapt_delta = 0.99))

Posterior <- posterior_samples(brm.student_stimuli, pars = c(
  "prior_Intercept",
  "b_Intercept",
  "prior_sd_study_ID",
  "sd_study_ID__Intercept",
  "prior_sd_study_ID:expt_condition",
  "sd_study_ID:expt_condition__Intercept",
  "sigma",
  "prior_nu",
  "nu"
))

plot(conditional_effects(brm.student_stimuli))
summary(brm.student_stimuli)

#Model 5:
brm.student_interaction <- 
  brm_multiple(data = MA_data_imp, family = student,
               hedge_g|se(se_hedge_g) ~ 1 + mean_age_1*test_lang + (1|study_ID/expt_condition),
               prior = priors2,
               sample_prior = T,
               file = "brm.student_interaction",
               iter = 20000, 
               warmup = 2000,
               chains = 2,
               cores = 2,
               control = list(adapt_delta = 0.99))

Posterior <- posterior_samples(brm.student_interaction, pars = c(
  "prior_Intercept",
  "b_Intercept",
  "prior_sd_study_ID",
  "sd_study_ID__Intercept",
  "prior_sd_study_ID:expt_condition",
  "sd_study_ID:expt_condition__Intercept",
  "sigma",
  "prior_nu",
  "nu"
))

pp_check(brm.student_interaction)
c_eff <- conditional_effects(brm.student_interaction, effects = 'mean_age_1:test_lang', spaghetti = T, nsamples = 150)
c_eff_plot <- plot(c_eff, mean = FALSE, points = T, point_args = c(alpha = 1, size = 3, show.legend = FALSE), spaghetti_args = c(alpha = 0.00001, size = 0.1), plot = FALSE)[[1]] +
  xlab("Mean Age in Days") +
  ylab("Hedges' g") +
  xlim(c(0, 400)) +
  ylim(c(-2, 3)) +
  ggtitle('Plot of Interaction between Mean Age & Test Language') +
  theme_bw()

c_eff_plot + scale_colour_manual(name = 'Test Language', 
                                 labels = c('Native', "Non-Native"), 
                                 values = c("native" = "grey32", "non-native" = "coral2")) +
  guides(color = guide_legend(override.aes = list(alpha = 0.5, size = 1))) +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=20))

summary(brm.student_interaction)

#Model 6:
brm.student_interaction_stimuli_age <- 
  brm_multiple(data = MA_data_imp, family = student,
               hedge_g|se(se_hedge_g) ~ 1 + mean_age_1*stimuli_type + (1|study_ID/expt_condition),
               prior = priors2,
               sample_prior = T,
               file = "brm.student_interaction_comp_age",
               iter = 20000, 
               warmup = 2000,
               chains = 2,
               cores = 2,
               control = list(adapt_delta = 0.99))

Posterior <- posterior_samples(brm.student_interaction_stimuli_age, pars = c(
  "prior_Intercept",
  "b_Intercept",
  "prior_sd_study_ID",
  "sd_study_ID__Intercept",
  "prior_sd_study_ID:expt_condition",
  "sd_study_ID:expt_condition__Intercept",
  "sigma",
  "prior_nu",
  "nu"
))

plot(conditional_effects(brm.student_interaction_stimuli_age))
summary(brm.student_interaction_stimuli_age)

#Check fits and model comparison:
brm.student_baseline <- add_criterion(brm.student_baseline, criterion="loo")
brm.student_age <- add_criterion(brm.student_age, criterion = "loo")
brm.student_age_mo <- add_criterion(brm.student_age_mo, criterion = "loo")
brm.student_lang <- add_criterion(brm.student_lang, criterion="loo")
brm.student_stimuli <- add_criterion(brm.student_stimuli, criterion="loo")
brm.student_interaction <- add_criterion(brm.student_interaction, criterion="loo")
brm.student_interaction_stimuli_age <- add_criterion(brm.student_interaction_stimuli_age, criterion="loo")

#loo model of baseline:
loo_model <- loo(brm.student_baseline)
#shows evidence of one outlier, i.e. pareto k-value is above 0.7:
plot(loo_model, label_points = T)
#algorithm for updating a loo object when Pareto k estimates are large.
loo_model <- loo_moment_match(brm.student_baseline, loo=loo_model)
plot(loo_model, label_points = T)
pareto_k_influence_values(loo_model)

#model comparison:
loo_model_weights(brm.student_baseline, brm.student_age)
loo_model_weights(brm.student_baseline, brm.student_lang)
loo_model_weights(brm.student_baseline, brm.student_stimuli)
loo_model_weights(brm.student_baseline, brm.student_interaction)
loo_model_weights(brm.student_baseline, brm.student_interaction_stimuli_age)

#Create forest plot:

#extract data from the baseline model:
study.draws <- spread_draws(brm.student_baseline, r_study_ID[study_ID,], b_Intercept) %>% 
  mutate(b_Intercept = r_study_ID + b_Intercept)

pooled.effect.draws <- spread_draws(brm.student_baseline, b_Intercept) %>% 
  mutate(study_ID = "Pooled Effect")

forest.data <- bind_rows(study.draws, pooled.effect.draws) %>% 
  ungroup() %>%
  mutate(study_ID = str_replace_all(study_ID, "[.]", " ")) %>% 
  mutate(study_ID = reorder(study_ID, b_Intercept))

forest.data.summary <- group_by(forest.data, study_ID) %>% 
  mean_qi(b_Intercept)

#plot the results:
ggplot(aes(b_Intercept, relevel(study_ID, "Pooled Effect", after = Inf)), 
       data = forest.data) +
  geom_vline(xintercept = 0, color = "black", size = 0.3, linetype="dotted") +
  geom_density_ridges(fill = "lightsteelblue4", rel_min_height = 0.03, col = NA, scale = 0.9,
                      alpha = 0.8, linetype = "dotted", size = 0.5) +
  geom_point(data=forest.data.summary, color='black', shape=18, size=2) +
  geom_errorbarh(data=forest.data.summary, aes(xmax = .upper, xmin = .lower, height = 0.1)) +
  geom_text(data = mutate_if(forest.data.summary, is.numeric, round, 3),
            aes(label = glue("{b_Intercept} [{.lower}, {.upper}]"), x = Inf), hjust = "inward", size = 3) +
  labs(x = "placeholder", 
       y = element_blank()) +
  scale_x_continuous("Hedges' g", limits = c(-0.75, 2.2)) +
  ggtitle("Forest Plot of Estimated Effect Sizes for Meta-Analytic Studies") +
  theme_classic()

#Examine publication bias: 

#create average of 100 imputations: use only for significance funnel plot & publication bias analysis
MA_data_imp_for_average <- mice(MA_data, meth=meth, post=post, print=FALSE, m=100, maxit=25, seed = 7)
all_data_imputations <- as_tibble(MA_data_imp_for_average$imp$se_hedge_g)
#mean across 100 data imputations:
mean_mi <- rowMeans(all_data_imputations)
mean_mi <- as_tibble(mean_mi)
#combine original with imputed sd values:
original_values <- na.omit(MA_data$se_hedge_g)

MA_data_average_imp <- MA_data %>%
  mutate(se_hedge_g = c(original_values, mean_mi$value)) %>%
  relocate(se_hedge_g, .before = n_1)

MA_data_average_imp <- MA_data_average_imp %>%
  mutate(vi = ((se_hedge_g )^2)) %>%
  relocate(vi, .before = se_hedge_g)

#calculate severity of publication bias needed to "explain away" the results:
svalue <- svalue( yi = MA_data_average_imp$hedge_g,
                  vi = MA_data_average_imp$vi,
                  q=0,
                  clustervar = MA_data_average_imp$study_ID,
                  model = "robust",
                  alpha.select = 0.05,
                  eta.grid.hi = 5,
                  favor.positive = TRUE,
                  CI.level = 0.95,
                  small = TRUE,
                  return.worst.meta = TRUE)

#make sensitivity plot, as in Mathur & VanderWeele (2020):
eta.list = as.list( c(150, rev( seq(1,100,1) ) ) )
res.list = lapply( eta.list, function(x) {
  cat("\n Working on eta = ", x)
  return( corrected_meta( yi = MA_data_average_imp$hedge_g,
                          vi = MA_data_average_imp$vi,
                          eta = x,
                          model = "robust",
                          favor.positive = TRUE,
                          clustervar = MA_data_average_imp$study_ID) )
}
)
# put results for each eta in a dataframe and plot:
res.df = as.data.frame( do.call( "rbind", res.list ) )
#plot the results:
Sensitivity_analysis <- ggplot( data = res.df, aes( x = eta, y = est ) ) + 
  geom_ribbon( data = res.df, aes( x = eta, ymin = lo, ymax = hi ), fill = "skyblue4" ) +
  geom_line( lwd = 1.5) +
  xlab( 'Publication Probability for Significant Studies') +
  ylab('Effect Size Esimate' ) + 
  geom_hline(yintercept = 0.0, linetype = "dotted", color = "orange", size = 0.7) +
  geom_hline(yintercept = svalue$meta.worst$b.r, linetype = "dashed", color = "red", size = 0.7) +
  scale_y_continuous(breaks=c(-0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) +
  scale_x_continuous(breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150)) +
  theme_classic()

Sensitivity_analysis + ggtitle("Sensitivity Analysis for Effect Size Estimate") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=20))

#create a significance funnel plot to examine publication bias:
sig_fun <- significance_funnel( yi = MA_data_average_imp$hedge_g,
                                vi = MA_data_average_imp$vi,
                                xmin = min(MA_data_average_imp$hedge_g),
                                xmax = max(MA_data_average_imp$hedge_g),
                                ymin = min(sqrt(MA_data_average_imp$se_hedge_g^2)) - 0.04,
                                ymax = max(sqrt(MA_data_average_imp$se_hedge_g^2)),
                                xlab = "Point Estimate of Effect Size",
                                ylab = "Sqaured Standard Error of Effect Size",
                                favor.positive = T,
                                alpha.select = 0.05)
sig_fun +
  ggtitle("Significance Funnel Plot of Meta-Analytic Studies") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=20)) +
  geom_point(aes(x=0.349, y=min(sqrt(MA_data_average_imp$vi)) - 0.05), size = 8, shape = 18, colour="black", show.legend = F) +
  geom_point(aes(x=0.237, y=min(sqrt(MA_data_average_imp$vi)) - 0.05), size = 8, shape = 18, colour="lightsteelblue4", show.legend = F) +
  geom_hline(yintercept = min(sqrt(MA_data_average_imp$vi)) - 0.05, linetype = "dashed")

#power calculations:
pwr.t.test(n = , d = svalue$meta.worst$b.r, sig.level = 0.05, power = 0.80, type = "one.sample")


#plot of posterior samples for pooled effect size:
post.samples <- posterior_samples(brm.student_baseline, c("^b", "^sd"))
names(post.samples) <- c("smd", "tau")
mean <- mean(post.samples$smd)

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


## Plot for the unique ID
Preds <- posterior_linpred(brm.student_baseline)
MA_data$unique <- NA
for (i in unique(MA_data$study_ID)){
  MA_data$unique[MA_data$study_ID==i] <- seq(nrow(MA_data[MA_data$study_ID==i,]))
}

MA_data <- MA_data %>% mutate(
  ID = paste0(study_ID,"_",unique)
)

d_pos <- data.frame(ID = NA, estimate = NA)

for (i in seq(length(unique(MA_data$ID)))){
  temp <- data.frame(ID = MA_data$ID[i], estimate = Preds[,i])
  d_pos <- rbind(d_pos, temp)
}

d_pos <- subset(d_pos, !is.na(ID))

pooledEffect <- spread_draws(brm.student_baseline, b_Intercept) %>% 
  mutate(ID = "Pooled Effect") %>%
  rename(estimate = b_Intercept) %>%
  select(ID, estimate)

d_pos <- rbind(d_pos, pooledEffect)

d_pos <- d_pos %>% mutate(
  ID = as.factor(ID),
  ID = relevel(ID, "Pooled Effect", after = Inf))

forest.data.summary <- group_by(d_pos, ID) %>% 
  mean_qi(estimate)



ggplot(aes(estimate, ID), 
       data = d_pos) +
  geom_vline(xintercept = 0, color = "black", size = 0.3, linetype="dotted") +
  geom_density_ridges(fill = "lightsteelblue4", rel_min_height = 0.03, col = NA, scale = 0.9,
                      alpha = 0.8, linetype = "dotted", size = 0.5) +
  geom_point(data=forest.data.summary, color='black', shape=18, size=2) +
  geom_errorbarh(data=forest.data.summary, aes(xmax = .upper, xmin = .lower, height = 0.1)) +
  geom_text(data = mutate_if(forest.data.summary, is.numeric, round, 3),
            aes(label = glue("{estimate} [{.lower}, {.upper}]"), x = Inf), hjust = "inward", size = 3) +
  labs(x = "placeholder", 
       y = element_blank()) +
  scale_x_continuous("Hedges' g", limits = c(-0.75, 2.2)) +
  ggtitle("Forest Plot of Estimated Effect Sizes for Meta-Analytic Studies") +
  theme_classic()
