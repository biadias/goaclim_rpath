library(MASS)

# Parameters for drawing distributions
  mu_1 <- 15  # Mean of log(sample type 1)
  mu_2 <- 15  # Mean of log(sample type 2) 

  cv_1   <-  1.0  # cv of sample type 1 (cv NOT in logspace)
  cv_2   <-  1.0  # cv of sample type 2 (cv NOT in logspace)
  cv_cov <- 0 #-0.2  # 0 for independent vars

# Specifying covar as a "cv", though cv is not an accurate way to describe it,
# scaling it this way helps scale it with respect to cv_1 and cv_2.
# Each of the covarying variables could have a different mean therefore
# different cv for the same actual covariance, so the number input below ends up 
# "between" sqrt(covariance)/mean1 and sqrt(covariance)/mean2. 
#  
# Note you will get an error in the cov matrix is not positive definite,
# to be positive definite, ensure that (sigma2_cov)^2 < (sigma2_1 * sigma2_2)

# Make a covariance matrix out of the cvs, using the variance and covariance 
# for generating log(samples) ~ Normal(mu, sigma2)) with specified covariance.
  sigma2_1   <- log(cv_1^2 + 1)
  sigma2_2   <- log(cv_2^2 + 1)
  sigma2_cov <- sign(cv_cov)*log(cv_cov^2 + 1)
  sigma_mat  <- diag(2)*c(sigma2_1, sigma2_2)
  sigma_mat[2,1] <- sigma_mat[1,2] <- sigma2_cov

# Number of tests to run
  Nsurveys  <- 1000
  
  
# Stations in a survey   
  Nstations <- 300 # trying similar ballpark to BTS

# draw the log(samples)
  log_samples <- mvrnorm(n=Nstations, mu=c(mu_1,mu_2), Sigma=sigma_mat, empirical=F)
# exponentiate to get to "actually measured" values, and also the sum of the 
# two variables.  Summing of the two variables can't be done directly in logspace.
  samples <- exp(log_samples)
  samples_stationsum <- rowSums(samples)

# now back-transform to do comparisons
  log_sumsamples <- log(samples_stationsum)

# collect statistics  
  mean_logsample1 <- mean(log_samples[,1])
  sd_logsample1   <- sd(log_samples[,1])
  se_logsample1   <- sd_logsample1/sqrt(Nstations)

  mean_logsample2 <- mean(log_samples[,2])
  sd_logsample2   <- sd(log_samples[,2])
  se_logsample2   <- sd_logsample2/sqrt(Nstations)
  
  mean_logsamplesum <- mean(log_sumsamples)
  sd_logsamplesum   <- sd(log_sumsamples)
  se_logsamplesum   <- sd_logsamplesum/sqrt(Nstations)

  mean_oneplustwo   <- mean_logsample1 + mean_logsample2
  # This is wrong
    se_oneplustwo     <- sqrt(se_logsample1^2 + se_logsample2^2)/sqrt()
  
  
  library(ggplot2)
  ggplot(data.frame(samples), aes(x=X1, y=X2) ) +
    geom_bin2d(bins = 70) +
    scale_fill_continuous(type = "viridis") +
    theme_bw()



  
