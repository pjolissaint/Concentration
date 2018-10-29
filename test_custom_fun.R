#============== SCRIPT FOR TESTING OF CUSTOM FUNCTIONS =================#


#========== Custom functions =========#
source("/home/pnj/R/R projects/Concentration/concentration_helper_fcts.R")
# Error function
ERFun <- function(x){
  2*pnorm(sqrt(2)*x) - 1
  }


#========== GINI Coefficients ==========#

## Test on Gini coeff
  #Tests
  test <- c()
  Gini_coeff_fun(test)
  HH_index_fun(test)
  
  test <- c(10, rep(0, 1e4))
  Gini_coeff_fun(test)
  HH_index_fun(test)
  test <- c(10, rep(0, 4))
  Gini_coeff_fun(test)
  HH_index_fun(test)
  test <- c(1000, rep(0, 4))
  Gini_coeff_fun(test)
  HH_index_fun(test)
  
  #
  k <- 17
  n <- 33
  test <- c(rep(1/k, k), rep(0, n-k))
  test <- c(rep(1, k), rep(0, n-k))
  Gini_coeff_fun(test)
  HH_index_fun(test)
  distr <- discr_distr_fun(test)
  lorenz <- discr_Lorenz_curve_fun(probas = distr$Prob, y = distr$Y, plot_curve = TRUE)
  
  
  test <- sample(x = c(10, 30, 50, 100, 200), size = 50, replace = TRUE)
  Gini_coeff_fun(test)
  HH_index_fun(test)
  distr <- discr_distr_fun(test)
  lorenz <- discr_Lorenz_curve_fun(probas = distr$Prob, y = distr$Y, plot_curve = TRUE)
  par(new = T)
  lines(x = seq(from = 0, to = 1, by = 1/length(test)), y = c(0, cumsum(sort(x = test, decreasing = FALSE))/sum(test)), col = "green")
  
  #unif distr
  a <- 7
  b <- 13
  test <- runif(n = 1e5, min = a, max = b)
  Gini_coeff_fun(test)
  # Unif - Theoretical
  (b-a)/(3*(b+a))
  
  lorenz <- discr_Lorenz_curve_fun(probas = distr$Prob, y = distr$Y, plot_curve = TRUE)
  
  
  # exponential distribution
  exp_mu <- 12
  test <- rexp(n = 1e5, rate = 1/exp_mu)
  Gini_coeff_fun(test) # should be 0.5 for any rate
  
  # Log normal
  mean_log <- 0
  sd_log <- 2
  test <- rlnorm(n = 1e6, meanlog = mean_log, sdlog = sd_log)
  Gini_coeff_fun(test)
  ERFun(sd_log/2) 
