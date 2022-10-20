library(Matrix)
library(tidyverse)
library(brms)



# Model comparison in Bayesian models

M_5 <- brm(weight ~ height, data = weight_df)
waic_M_3 <- waic(M_3)
waic_M_5 <- waic(M_5)
loo_compare(waic_M_3, waic_M_5)

loo_M_3 <- loo(M_3)
loo_M_5 <- loo(M_5)
loo_compare(loo_M_3, loo_M_5)

# emmeans
library(emmeans)
ref_grid(M_3)
emmeans(ref_grid(M_3), specs = ~ gender)
emmeans(ref_grid(M_3), specs = pairwise ~ gender)


# Non-homogeneity of variance ---------------------------------------------

set.seed(10101)
n <- 250
data_df2 <- tibble(A = rnorm(n, mean = 1, sd = 1),
                   B = rnorm(n, mean = 0.25, sd = 2)
) %>% pivot_longer(cols = everything(), names_to = 'x', values_to = 'y')

ggplot(data_df2, aes(x,y)) + geom_boxplot()

M_6 <- brm(y ~ x, data = data_df2)
M_7 <- brm(
  bf(y ~ x, sigma ~ x),
  data = data_df2
)
M_7


# Robust regression -------------------------------------------------------

set.seed(10101)
n <- 25
data_df3 <- tibble(x = rnorm(n),
                   y = 5 + 0.5 * x + rnorm(n, sd = 0.1))


data_df4 <- data_df3
data_df4[12,2] <- 7.5

ggplot(data_df3, aes(x, y)) + geom_point()
ggplot(data_df4, aes(x, y)) + geom_point()

# frequentist analysis 
M_freq_3 <- lm(y ~ x, data = data_df3)
M_freq_4 <- lm(y ~ x, data = data_df4)

summary(M_freq_3)
summary(M_freq_4)

# Bayesian normal linear model on outlier data
M_8 <- brm(y ~ x, data = data_df4)
# Bayesian non-normal linear model on outlier data
M_9 <- brm(y ~ x, 
           family = student(),
           data = data_df4)

mcmc_plot(M_8, type = 'areas') 
mcmc_plot(M_9, type = 'areas', variable = c('b_Intercept', 'b_x', 'sigma'))
M_9


# Bayesian binary logistic regression -------------------------------------

weight_df1 <- mutate(weight_df, dweight = weight > median(weight))

M_10 <- brm(dweight ~ height + age + gender,
            data = weight_df1,
            family = bernoulli())

M_freq_5 <- glm(dweight ~ height+ age + gender,
                data = weight_df1,
                family = binomial())

M_10
fixef(M_10)
summary(M_freq_5)$coefficients


# Bayesian generalized linear models for counts ---------------------------

# frequentst Poisson regression
M_freq_6 <- glm(publications ~ prestige,
                data = biochemists_df,
                family = poisson())

# Poisson regression
M_11 <- brm(publications ~ prestige,
            data = biochemists_df,
            family = poisson())

# negative binomial regression
M_12 <- brm(publications ~ prestige,
            data = biochemists_df,
            family = negbinomial())

# model comparison
waic_M_11 <- waic(M_11)
waic_M_12 <- waic(M_12)
loo_compare(waic_M_11, waic_M_12)
waic_11_12 <- waic(M_11, M_12)
-2 * waic_11_12$diffs


# predicted mean of the the distribution for any given value of `prestige`
prestige_df <- tibble(prestige = seq(5))
# get mean of posterior over predicted log mean for each value of `prestige`
posterior_linpred(M_11, newdata = prestige_df) %>% apply(2, mean)
# get percentiles (0, 25, 50, 75, 100) for posterior over log mean for each value of `prestige`
posterior_linpred(M_11, newdata = prestige_df) %>% apply(2, quantile)

# get percentiles (0, 25, 50, 75, 100) for posterior over mean for each value of `prestige`
posterior_linpred(M_11, newdata = prestige_df, transform = TRUE) %>% apply(2, quantile)

# get 95% posterior interval of posterior over mean for each value of `prestige`
hpd <- function(x)quantile(x, probs = c(0.025, 0.975))
posterior_linpred(M_11, newdata = prestige_df, transform = TRUE) %>% apply(2, hpd)

# And for the Negative Binomial regression  .....
posterior_linpred(M_12, newdata = prestige_df, transform = TRUE) %>% apply(2, hpd)




# Zeroinflated count data -------------------------------------------------

cig_mean <- mean(smoking_df$cigs)
# what is probability of zero in Poisson distribution with mean of `cig_mean`
dpois(x = 0, lambda = cig_mean)


# Poisson regression
M_13 <- brm(cigs ~ educ, 
            data = smoking_df,
            family = poisson())

# Neg binomial regression
M_14 <- brm(cigs ~ educ, 
            data = smoking_df,
            family = negbinomial())


# Zero inflated Poisson regression
M_15 <- brm(cigs ~ educ, 
            data = smoking_df,
            family = zero_inflated_poisson())

# Zero inflated Neg binomial regression
M_16 <- brm(cigs ~ educ, 
            data = smoking_df,
            family = zero_inflated_negbinomial())


# Zero inflated Poisson regression ++ 
M_17 <- brm(bf(cigs ~ educ, zi ~ educ),
            data = smoking_df,
            family = zero_inflated_poisson())

# Zero inflated Neg binomial regression ++ 
M_18 <- brm(bf(cigs ~ educ, zi ~ educ), 
            data = smoking_df,
            family = zero_inflated_negbinomial())

zi_waic <- waic(M_13, M_14, M_15, M_16, M_17, M_18)
-2 * zi_waic$diffs


# Mixed effects models -----------------------------------------------------
library(lme4)
ggplot(sleepstudy,
       aes(x = Days, y = Reaction, colour = Subject)
) + geom_point() + 
  stat_smooth(method = 'lm', se = F) +
  facet_wrap(~Subject)

# frequentist linear mixed effects model
M_freq_7 <- lmer(Reaction ~ Days + (Days|Subject), data = sleepstudy)
summary(M_freq_7)

# frequentist linear mixed effects model
M_19 <- brm(Reaction ~ Days + (Days|Subject), data = sleepstudy)
