#'manual_logistic_regression
#'
#'This is an function named 'manual_logistic_regression' which xxx
#'
#'@param x input value
#'
#'@return the square of x
#'
#'@examples
#'square(3)
#'
#'@export
#'


# Newton-Raphson Method
manual_logistic_regression <- function(x, y, threshold = 1e-10, max_iter = 100) {

  #A function to return p, given X and beta
  #We'll need this function in the iterative section
  calc_p = function(x,beta) {
    beta = as.vector(beta)
    p = exp(x%*%beta) / (1 + exp(x%*%beta))
    return(p)
  }

  #### setup bit ####
  #make sure x is in matrix form
  x = as.matrix(x)

  #bias in order to calculate the intercept
  x0 = rep(1,30)
  x = cbind(x0,x)

  #initial guess for beta
  beta = rep(0,ncol(x))

  #initial value bigger than threshold so that we can enter our while loop
  diff = 10000

  #counter to ensure we're not stuck in an infinite loop
  iter_count = 0

  #### iterative bit ####
  while(diff > threshold) { #tests for convergence

    #calculate probabilities using current estimate of beta
    p = as.vector(calc_p(x,beta))

    #calculate matrix of weights W
    w =  diag(p*(1-p))

    #calculate the change in beta
    beta_change = solve(t(x)%*%w%*%x) %*% t(x)%*%(y-p)

    #update beta
    beta = beta + beta_change

    #calculate how much we changed beta by in this iteration
    #if this is less than threshold, we'll break the while loop
    diff = sum(beta_change^2)

    #see if we've hit the maximum number of iterations
    iter_count = iter_count + 1
    if(iter_count > max_iter) {
      break
    }
  }

  # parameter coefficients
  coef = c("(Intercept)" = beta[1], x = beta[-1])

  # number of iteration
  fisher_scoring_iterations = iter_count

  # variance-covariance matrix
  covariance_matrix = solve(t(x)%*%w%*%x)

  # Dispersion parameter for binomial family taken to be 1
  dispersion = 1

  # standard error
  se = sqrt(diag(covariance_matrix)*dispersion)

  # test statistics
  t = beta/se
  tstat = c("(Intercept)" = t[1], x = t[-1])

  # p-value
  pvalue = 2*(1-pnorm(abs(t)))
  p_value = c("(Intercept)" = pvalue[1], x = pvalue[-1])

  # Deviance Residuals
  # The calculate predicted probabilities
  pred_p = exp(x%*%beta) / (1 + exp(x%*%beta))
  # Apply formula for deviance residuals
  dev_res = sign(y-pred_p) * ifelse(y==1,sqrt(-2*log(pred_p)),sqrt(-2*log(1-pred_p)))

  # summary of deviance residuals
  sum_dev_res = summary(dev_res)

  # Residual deviance
  res_dev = -2*sum(y*log(pred_p) + (1 - y)*log(1 - pred_p))

  # Residual deviance degrees of freedom
  df_residual = dim(x)[1] - dim(x)[2]

  # Null deviance
  null_dev = -2*sum(y*log(mean(y)) + (1 - y)*log(1 - mean(y)))

  # Null deviance degrees of freedom
  df_null = dim(x)[1] - 1

  # Akaike Information Criterion
  AIC = 2*ncol(x) + res_dev

  z = list(coefficients = coef,
           covariance_matrix = covariance_matrix,
           Dispersion = dispersion,
           standard_error = se,
           t_statistics = tstat,
           p_value = p_value,
           deviance_residuals = dev_res,
           quantile_deviance_residuals = sum_dev_res,
           Null_deviance_degrees_of_freedom = df_null,
           Residual_deviance_degrees_of_freedom = df_residual,
           Null_deviance = null_dev,
           Residual_deviance = res_dev,
           AIC = AIC,
           number_iterations = fisher_scoring_iterations
  )

  return(z)

}

set.seed(2022)
#simulate data
#independent variables
x1 = rnorm(30,3,2) + 0.1*c(1:30)
x2 = rbinom(30, 1,0.3)
x3 = rpois(n = 30, lambda = 4)
x3[16:30] = x3[16:30] - rpois(n = 15, lambda = 2)

#dependent variable
y = c(rbinom(5, 1,0.1),rbinom(10, 1,0.25),rbinom(10, 1,0.75),rbinom(5, 1,0.9))

x = cbind(x1,x2,x3)

manual_logistic_regression(x,y)
summary(glm(y~x1+x2+x3, family = "binomial"))

#vcov(glm(y~x1+x2+x3, family = "binomial"))

# McFadden's pseudo R^2
#  r2Log = 1 - dev_res / dev_null
#           R_squared = r2Log,

# p-value for the chi^2
#  p_value_chi = 1 - pchisq(res_dev, df_residual)
