library(Matrix)
library(tidyverse)
library(brms)


# Generate some data ------------------------------------------------------

set.seed(10101) # Omit or change this if you like

N <- 25

x_1 <- rnorm(N)
x_2 <- rnorm(N)

beta_0 <- 1.25
beta_1 <- 1.75
beta_2 <- 2.25

mu <- beta_0 + beta_1 * x_1 + beta_2 * x_2

y <- mu + rnorm(N, mean=0, sd=1.75)

data_df1 <- tibble(x_1, x_2, y)


# Frequentist linear model ------------------------------------------------

M_freq_1 <- lm(y ~ x_1 + x_2, data = data_df1)
# look at coefs
coef(M_freq_1)
# look at residual standard dev
sigma(M_freq_1)
# get summary 
summary(M_freq_1)

# Bayesian linear model
M_1 <- brm(y ~ x_1 + x_2, data = data_df1)
M_1

# default plots
plot(M_1)

# histograms of posteriors
mcmc_plot(M_1, type = 'hist', binwidth = 0.05)
# plot posterior over intercept only
mcmc_plot(M_1, variable = 'b_Intercept', type = 'hist', binwidth = 0.05)
# plot posterior over slope terms only
mcmc_plot(M_1, variable = 'b_x.*', regex = TRUE, type = 'hist', binwidth = 0.05)

# area plot
mcmc_plot(M_1, type = 'areas')
mcmc_plot(M_1, type = 'areas_ridges', prob = 0.5)

# look at the priors
prior_summary(M_1)
# see priors in a model before it is run
get_prior(y ~ x_1 + x_2, data = data_df1)


# changing defaults -------------------------------------------------------

M_2 <- brm(y ~ x_1 + x_2, 
           iter = 2500,
           warmup = 500,
           chains = 6,
           cores = 4,
           prior = set_prior('normal(0, 100)'),
           data = data_df1)

# confirm priors are as expected
prior_summary(M_2)

# look at coefficients only
fixef(M_2)

# compare to the M_1 case, with different priors
fixef(M_1)



# Read in all the data sets -----------------------------------------------

biochemists_df <- read_csv('https://raw.githubusercontent.com/mark-andrews/bada02/master/data/biochemist.csv')
smoking_df <- read_csv('https://raw.githubusercontent.com/mark-andrews/bada02/master/data/smoking.csv')
classroom_df <- read_csv("https://raw.githubusercontent.com/mark-andrews/bada02/master/data/classroom.csv")
mathach_df <- read_csv("https://raw.githubusercontent.com/mark-andrews/bada02/master/data/mathachieve.csv")
mathach_school_df <- read_csv <- read_csv("https://raw.githubusercontent.com/mark-andrews/bada02/master/data/mathachieveschool.csv")
science_df <- read_csv('https://raw.githubusercontent.com/mark-andrews/bada02/master/data/science.csv')
owls <- read_csv('https://raw.githubusercontent.com/mark-andrews/bada02/master/data/owls.csv')
weight_df <- read_csv("https://raw.githubusercontent.com/mark-andrews/bada02/main/data/weight.csv")

# linear regression predicting weight from gender and height and age
M_freq_2 <- lm(weight ~ height + age + gender, data = weight_df)

# Bayesian version using default priors
M_3 <- brm(weight ~ height + age + gender, data = weight_df)

prior_summary(M_3)

new_priors <- c(
  set_prior(class = 'b', coef = 'age', prior = 'normal(0, 10)'),
  set_prior(class = 'b', coef = 'gendermale', prior = 'normal(0, 10)'),
  set_prior(class = 'b', coef = 'height', prior = 'normal(0, 10)'),
  set_prior(class = 'Intercept', prior = 'normal(0, 100)'),
  set_prior(class = 'sigma', prior = 'student_t(1, 0, 25)')
)

# Bayesian version using default priors
M_4 <- brm(weight ~ height + age + gender, 
           data = weight_df,
           prior = new_priors)

# look at coefficients
round(fixef(M_4), 3)
# compare with lm model
round(summary(M_freq_2)$coefficients, 3)
round(confint(M_freq_2), 3)

summary(M_freq_2)$r.sq
bayes_R2(M_4)
