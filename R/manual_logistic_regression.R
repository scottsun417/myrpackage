#'manual_logistic_regression
#'
#'This is an function named 'manual_logistic_regression' which implements a logistic regression
#'used for modelling categorical data and returns a list containing sufficient model summary information.
#'This function assuming the response variable Y falls into exactly two categories, Y = 1 and
#'the other by Y = 0.
#'
#'@param x input value, predictors preferred in a matrix form
#'@param y input value, response vector consists only of 0 and 1
#'@param threshold input value, a numeric value checking if the logistic regression algorithm is
#'converging or not (the amount of beta changes should get smaller and smaller with each iteration,
#'eventually falling below the threshold), default threshold = 1e-10
#'@param max_iter input value, the maximum number of iterations, default max_iter = 100
#'
#'@return A list containing model summary information, including:
#'@return coefficients: estimated parameters for intercept and each variable (predictor),
#'according to the response
#'@return covariance_matrix: used to derive the standard errors and confidence intervals of the
#'fitted model's coefficient estimates
#'@return Dispersion: dispersion parameter for binomial family taken to be 1
#'@return standard_error: standard error for intercept and each variable (predictor),
#'according to the response
#'@return t_statistics: test statistics for intercept and each variable (predictor),
#'according to the response
#'@return p_value: P-value for intercept and each variable (predictor), according to the response
#'@return deviance_residuals: measuring how much probabilities estimated from fitted model differ
#'from the observed proportions of successes
#'@return quantile_deviance_residuals: a five-number summary of the residuals, including Minimum,
#'1st Quantile, Median, 3rd Quantile and Maximum, respectively
#'@return Null_deviance_degrees_of_freedom: the degrees of freedom of Null deviance
#'@return Residual_deviance_degrees_of_freedom: the degrees of freedom of Residual deviance
#'@return Null_deviance: measuring how well the response variable can be predicted by a model with
#'only an intercept term
#'@return Residual_deviance: measuring how well the response variable is predicted by the model when the
#'predictors are included
#'@return AIC: Akaike Information Criterion
#'@return number_iterations: the number of iterations
#'
#'@importFrom stats fivenum
#'@importFrom stats pnorm
#'
#'@examples
#'set.seed(2022)
#'x1 = rnorm(30,3,2) + 0.1*c(1:30)
#'x2 = rbinom(30, 1,0.3)
#'x3 = rpois(n = 30, lambda = 4)
#'x3[16:30] = x3[16:30] - rpois(n = 15, lambda = 2)
#'x = cbind(x1,x2,x3)
#'y = c(rbinom(5, 1,0.1),rbinom(10, 1,0.25),rbinom(10, 1,0.75),rbinom(5, 1,0.9))
#'manual_glm <- manual_logistic_regression(x, y)
#'
#'@export
#'
# Newton-Raphson Method
manual_logistic_regression <- function(x, y, threshold = 1e-10, max_iter = 100) {

  # a function to return p, given x and beta
  # we will need this function in the iteration section
  calc_p = function(x,beta) {
    beta = as.vector(beta)
    p = exp(x%*%beta) / (1 + exp(x%*%beta))
    return(p)
  }

  #### Setup ####

  # make sure x in matrix form
  x = as.matrix(x)

  # bias in order to calculate the intercept
  x0 = rep(1,30)
  x = cbind(x0,x)

  # initial guess for beta
  beta = rep(0,ncol(x))

  # initial value bigger than threshold so that we can do the while loop
  diff = 10000

  # counter to ensure not in an infinite loop during iterative
  iter_count = 0

  #### Iteration ####

  # do when convergence (check converging or not)
  while(diff > threshold) {

    # calculate probabilities using current estimate of beta
    p = as.vector(calc_p(x,beta))

    # calculate matrix of weights w
    w =  diag(p*(1-p))

    # calculate the change in beta
    beta_change = solve(t(x)%*%w%*%x) %*% t(x)%*%(y-p)

    # update beta
    beta = beta + beta_change

    # calculate how much we changed beta by in this iteration
    # if this is less than threshold, we will break the while loop
    diff = sum(beta_change^2)

    # see if it reach the maximum number of iterations
    iter_count = iter_count + 1
    if(iter_count > max_iter) {
      return(cat("Not Converging"))
    }
  }

  #### Outcome ####

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
  # the calculate predicted probabilities
  pred_p = exp(x%*%beta) / (1 + exp(x%*%beta))
  # apply formula for deviance residuals
  dev_res = sign(y-pred_p) * ifelse(y==1,sqrt(-2*log(pred_p)),sqrt(-2*log(1-pred_p)))

  # summary of deviance residuals
  sum_dev_res = fivenum(dev_res)

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

  #### Return list ####

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
