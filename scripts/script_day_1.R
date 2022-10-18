library(priorexposure)

# likelihood function for the case of 
# 139 heads in 250 flips
bernoulli_likelihood(n = 250, m = 139)

# likelihood function for the case of 
# 14 heads in 25 flips
bernoulli_likelihood(n = 25, m = 14)

# Plot a beta distribution with alpha = 3, beta = 5
beta_plot(alpha = 3, beta = 5)

# Beta plot centred at 0.5, and quite narrow
beta_plot(alpha = 10, beta = 10)
beta_plot(alpha = 20, beta = 20)
beta_plot(alpha = 100, beta = 100)

# Beta plot centred at 0.5, and less narrow
beta_plot(alpha = 5, beta = 5)
beta_plot(alpha = 2, beta = 2)
beta_plot(alpha = 1.0, beta = 1.0)

# non-symmetrical beta plot
beta_plot(alpha = 5, beta = 10)

# beta plots with hyperpars less than 1.0
beta_plot(alpha = 0.5, beta = 0.5)

# if m = 139, and n = 250
# and we choose alpha = beta = 5 as prior
# then, the posterior distribution is this:
beta_plot(5 + 139, 5 + 250 - 139)
bernoulli_posterior_plot(n = 250, 
                         m = 139, 
                         alpha = 5, 
                         beta = 5)

# summarize the posterior distribution 
# that we just plotted
bernoulli_posterior_summary(n = 250,
                            m = 139,
                            alpha = 5,
                            beta = 5)

# get 95% HPD interval
get_beta_hpd(5 + 139, 5 + 111)

# what if we used alpha = beta = 1 as prior
# get 95% HPD interval
get_beta_hpd(1 + 139, 1 + 111)
