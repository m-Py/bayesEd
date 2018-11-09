
## Author: Martin Papenberg, inspired from a blog post by Jeff Rouder 
## (http://jeffrouder.blogspot.com/2016/01/what-priors-should-i-use-part-i.html) 

## This file contains functions to visualize Bayes factors (prior and
## marginal likelihoods), and computes t-test Bayes factors for the
## comparison of two group means. The code is meant for educational
## purposes only; to actually compute Bayes factors, I encourage the use
## of the `BayesFactor` package.

#' Visualize a prior on Cohen's d for a Bayesian t-test
#' 
#' @param alternative A function object. The default is a Cauchy prior
#'     with scaling parameter `sqrt(2) / 2` as is the default in package
#'     `BayesFactor` (Morey & Rouder, 2015). This argument can also be a
#'     scalar number, in which case it is assumed that the alternative
#'     is a point hypothesis on Cohen's d (with Cohen's d = `prior`).
#' @param null Boolean. Plot a line on x = 0 representing the null
#'     hypothesis?
#' @param from the left margin of the x-axis
#' @param to the right margin of the x-axis
#' @param xlab the label of the x-axis
#' @param frame.plot Should a frame be drawn?
#' @param col The color of the curve of the alternative hypothesis.
#' @param lwd The line width
#' @param ... additional parameters passed to function `curve`
#' 
#' @examples
#' ## Standard cauchy prior:
#' visualize_prior(lty = 2, lwd = 2)
#' ## Wide cauchy prior:
#' visualize_prior(function(x) dcauchy(x, scale = 1))
#' ## Normal prior with M = 0 and SD = 0.4
#' visualize_prior(function(x) dnorm(x, 0, 0.4))
#' ## Point prior with Cohen's d = 0.4
#' visualize_prior(0.4)
#'
#' @references
#'
#' Morey, R. D., & Rouder, J. N. (2015). BayesFactor: Computation of
#'     bayes factors for common designs. Retrieved from
#'     https://CRAN.R-project.org/package=BayesFactor
#'
#' Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D., & Iverson,
#'     G. (2009). Bayesian t tests for accepting and rejecting the null
#'     hypothesis. Psychonomic Bulletin & Review, 16(2), 225–237.
#'
#' @export
#' 
visualize_prior <- function(alternative = function(x) dcauchy(x, scale = sqrt(2) / 2),
                            null = TRUE, from = -4, to = 4,
                            xlab = "Cohen's d", frame.plot = FALSE,
                            col = "red", lwd = 3, ...) {
  if (class(alternative) != "function") {
    plot(c(alternative, alternative), c(0, 1), type = "l",
         las = 1, xlim = c(from, to), xlab = xlab, 
         frame.plot = frame.plot, yaxt = "n", col = "transparent", 
         ylab = "")
    abline(v = alternative, col = col, lwd = lwd)
  } else {
    curve(alternative, las = 1, from = from, to = to, xlab = xlab, 
          frame.plot = frame.plot, yaxt = "n", col = col, ylab = "", 
          lwd = lwd, ...)
  }
  if (null) {
    abline(v = 0, col = "#005A31", lwd = lwd)
    legend("topleft", lwd = lwd, legend = c("Null", "Alternative"), 
           cex = 1.2, col = c("#005A31", "red"), lty = c(2,2), bty = "n")
  }
}


#' Visualize hypotheses' predictions for a Bayesian two-group t-test
#' 
#' Visualizes the predictions of the null and an alternative hypothesis.
#' By specifying an observed effect (observed t-value), an illustration
#' of the ratio of the marginal likelihoods, i.e., an
#' illustration of the Bayes factor is displayed.
#' 
#' @param alternative A function object. The default is a Cauchy prior
#'     with scaling parameter `sqrt(2) / 2` as is the default in package
#'     `BayesFactor` (Morey & Rouder, 2015). This argument can also be a
#'     scalar number, in which case it is assumed that the alternative
#'     is a point hypothesis on Cohen's d (with Cohen's d = `prior`).
#'     If `alternative` is NULL, only the null hypothesis predictions
#'     are drawn.
#' @param n1 The sample size in group 1
#' @param n2 The sample size in group 2
#' @param from the left margin of the x-axis
#' @param to the right margin of the x-axis
#' @param lty The line type; lty[1] specifies the line type for the null
#'     hypothesis; lty[2] specifies the line type for the alternative
#'     hypothesis; lty[3] specifies the line type for the alternative
#'     hypothesis (only has an effect if an observed effect size is
#'     specified via argument `observed_t`)
#' @param col Specify the coloring of the plot. col[1] specifies the
#'     color for the null hypothesis; lty[2] specifies the color for the
#'     alternative hypothesis; lty[3] specifies the color for the line
#'     the illustrates the observed effect (only has an effect if an
#'     observed effect size is specified via argument `observed_t`);
#'     lty[4] specifies the color of the points at the intersection of
#'     the observed effect and the curves for the hypotheses (only has
#'     an effect if an observed effect size is specified via argument
#'     `observed_t`)
#' @param observed_t An observed t-value. Does not need to be specified
#'     and defaults to NULL. If this argument is passed, the marginal
#'     likelihoods of the null and the alternative hypothesis are drawn
#'     and a Bayes factor is displayed.
#' @param frame.plot Should a frame be drawn?
#' @param xlab The label of the x-axis
#' @param BFx The x coordinate where the Bayes factor is displayed in
#'     the plot. Only has an effect if an `observed_t` is passed.
#' @param BFy The y coordinate where the Bayes factor is displayed in
#'     the plot. Only has an effect if an `observed_t` is passed.
#' @param BF10 Boolean; defaults to TRUE. Should a Bayes factor be
#'     displayed quantifying the evidence in favor of the alternative
#'     hypothesis? The BF01 is displayed if `BF10` is `FALSE`.  Only has
#'     an effect if an `observed_t` is passed.
#' @param legend_placement A keyword passed to `legend` to indicate
#'     where the legend has to be placed. Defaults to "topleft".  Only
#'     has an effect if an `observed_t` is passed.
#' @param ... additional parameters passed to `curve`
#'
#' @references
#'
#' Morey, R. D., & Rouder, J. N. (2015). BayesFactor: Computation of
#'     bayes factors for common designs. Retrieved from
#'     https://CRAN.R-project.org/package=BayesFactor
#'
#' Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D., & Iverson,
#' G. (2009). Bayesian t tests for accepting and rejecting the null
#' hypothesis. Psychonomic Bulletin & Review, 16(2), 225–237.
#' 
#' @examples 
#' ## Using the default cauchy prior
#' visualize_predictions(n1 = 100, n2 = 100, observed_t = 3, lty = c(1, 1, 2))
#' ## Using a wide Cauchy prior
#' visualize_predictions(function(x) dcauchy(x, scale = 1), 
#' n1 = 100, n2 = 100, observed_t = 3, lty = c(1, 1, 2))
#' ## Using a normal prior (M = 0, SD = 0.4)
#' visualize_predictions(function(x) dnorm(x, 0, 0.4), 
#' n1 = 100, n2 = 100, observed_t = 3, lty = c(1, 1, 2))
#'
#' @export
#' 
visualize_predictions <- function(alternative = function(x) dcauchy(x, scale = sqrt(2) / 2),
                                  n1, n2, from = -6, to = 6,
                                  lty = c(1, 1, 2), lwd = 3,
                                  col = c("#005A31", "red", "black", "orange"),
                                  observed_t = NULL, frame.plot = FALSE,
                                  xlab = "t-value", BFx = from + 0.1,
                                  BFy = 0.1, BF10 = TRUE,
                                  legend_placement = "topleft", ...) {
  ## 1. Draw null predictions
  nullhyp <- function(x) dt(x, n1 + n2 - 2)
  curve(nullhyp, from = from , to = to, col = col[1], lty = lty[1], lwd = lwd[1],
        frame.plot = frame.plot, yaxt = "n", ylab = "", xlab = xlab, ...)
  maxdensity <- nullhyp(0)
  ## 2. Draw alternative hypothesis predictions
  if (!is.null(alternative)) {
    predictions_alt <- sample_likelihoods(alternative, n1, n2, from, to)
    points(predictions_alt$x, predictions_alt$y, type = "l", col = col[2],
           lty = lty[2], lwd = lwd)
    legend(legend_placement, lwd = 3, legend = c("Null", "Alternative"), 
           cex = 1.2, col = c("#005A31", "red"), lty = c(2,2), bty = "n")
    maxdensity <- max(maxdensity, predictions_alt$y)
  }
  ## 3. Illustrate the observed effect and the Bayes factor
  if (!is.null(observed_t)) {
    null_likelihood <- nullhyp(observed_t)
    lines(x = rep(observed_t, 2), y = c(0, maxdensity), 
          col = col[3], lwd = lwd, lty = lty[3])
    points(observed_t, null_likelihood, pch = 19, col = col[4], cex = 2)
    if (!is.null(alternative)) {
      alt_likelihood <- marginal_likelihood(observed_t, alternative, n1, n2)
      BF <- alt_likelihood / null_likelihood
      displayBF <- "BF10 = "
      if (BF10 == FALSE) {
        BF <- 1 / BF
        displayBF <- "BF01 = "
      }
      points(observed_t, alt_likelihood, pch = 19, col = col[4], cex = 2)
      max_likelihood <- max(null_likelihood, alt_likelihood)
      text(x = BFx, y = BFy, pos = 4, 
           labels = paste(displayBF, round(BF, 2)))
    }
  }
}

#' Generate marginal likelihoods in a given interval of t-values
#' 
#' @param alternative A function object. The default is a Cauchy prior
#'     with scaling parameter `sqrt(2) / 2` as is the default in package
#'     `BayesFactor` (Morey & Rouder, 2015). This argument can also be a
#'     scalar number, in which case it is assumed that the alternative
#'     is a point hypothesis on Cohen's d (with Cohen's d = `prior`).
#' @param n1 The sample size in group 1
#' @param n2 The sample size in group 2
#' @param from The lowest t-value
#' @param to The largest t-value
#' 
#' @return A data.frame of two columns: Column x contains the t-values, 
#'   column y contains the marginal likelihoods
#'
#' @examples
#' sample_likelihoods(n1 = 30, n2 = 30, from = -3, to = 3)
#' sample_likelihoods(function(x) dnorm(x, 0, 0.3), n1 = 30, n2 = 30)
#' 
#' @export
#' 
sample_likelihoods <- function(alternative = function(x) dcauchy(x, scale = sqrt(2) / 2), 
                               n1, n2, from = -6, to = 6) {
  observed_ts <- seq(from, to, .01) # hypothetical observed t-values
  ## Allow for a point null hypothesis
  if (class(alternative) != "function") {
    althyp <- function(x) dt(x, df = n1 + n2 - 2, ncp = ncp_from_d(alternative, n1, n2))
    predictions <- althyp(observed_ts)
    return(data.frame(x = observed_ts, y = predictions))
  }
  predictions <- vector(length = length(observed_ts)) # predict marginal likelihood of each t
  for (i in 1:length(predictions)) {
    predictions[i] <- marginal_likelihood(observed_ts[i], alternative, n1, n2)
  }
  return(data.frame(x = observed_ts, y = predictions))
}

#' Compute the marginal likelihood of an observed t-value
#' 
#' This function can be used to compute a Bayes factor as the ratio
#' of two marginal likelihoods.
#' 
#' @param observed_t The observed t-value.
#' @param prior A function object describing the a-priori plausibility
#'     of all effect sizes. The default is a Cauchy prior as in package
#'     BayesFactor (Morey & Rouder, 2015). Can also be a scalar, in
#'     which case it is assumed that the alternative is a point
#'     hypothesis on Cohen's d (with value `prior`).
#' @param n1 The sample size in group 1
#' @param n2 The sample size in group 2
#'
#' @return The marginal likelihood of the observed t-value given the
#'     prior on Cohen's d.
#'
#' @examples
#' ## Compute Bayes factor for standard Cauchy prior:
#' n1 <- 100
#' n2 <- 100
#' sample1 <- rnorm(n1, 0.2, 1)
#' sample2 <- rnorm(n1, 0.2, 1)
#' tvalue  <- t.test(sample1, sample2)$statistic
#' ml_alt  <- marginal_likelihood(tvalue, n1 = n1, n2 = n2)
#' bf10    <- ml_alt / dt(tvalue, n1 + n2 - 2)
#' # Compare:
#' BayesFactor::ttestBF(sample1, sample2)
#' 
#' 
#' ## Normal prior - N(0, 0.3)
#' marginal_likelihood(3, function(x) dnorm(x, 0, 0.3), n1 = 30, n2 = 30)
#'
#' @export
#' 
marginal_likelihood <- function(observed_t, 
                                prior = function(x) dcauchy(x, scale = sqrt(2) / 2), 
                                n1, n2) {
  ## allow for point alternative:
  if (class(prior) != "function") {
    althyp <- function(x) dt(x, df = n1 + n2 - 2, ncp = ncp_from_d(prior, n1, n2))
    return(althyp(observed_t))
  }
  integrate(weighted_likelihood, lower = -Inf, upper = Inf, 
            observed_t = observed_t, prior = prior, n1 = n1, n2 = n2)$value
}

#' Compute the weighted likelihood of an observed t-value
#' 
#' For an assumed true effect size (Cohen's d), compute the likelihood
#' of an observed t-value weighted by the a priori plausbility of this
#' effect size. This function is primarily used to compute the marginal 
#' likelihood of the effect size given a continuous hypothesis (e.g., 
#' a Cauchy prior) by integrating in function `marginal_likelihood`
#' below.
#' 
#' @param true_d A one-element vector, representing an assumed true 
#'   effect size d, conditioned on which the likelihood is computed.
#' @param observed_t The observed t-value.
#' @param prior A function object describing the a-priori plausibility
#'   of all effect sizes. The default is a cauchy prior
#'   as in package BayesFactor (Morey & Rouder, 2015)
#' @param n1 The sample size in group 1
#' @param n2 The sample size in group 2
#'
#' @return The weighted likelihood of the observed t-value.
#'
#' @export
#' 
weighted_likelihood <- function(true_d, observed_t, 
                                prior = function(x) dcauchy(x, scale = sqrt(2) / 2), 
                                n1, n2) {
  dfr <- n1 + n2 - 2 ## degrees of freedom for two-sample t-test
  ncp = ncp_from_d(true_d, n1, n2)
  return(dt(observed_t, dfr, ncp) * prior(true_d))
}

                     
#' Compute non-centrality parameter for the non-central t distribution
#' 
#' @param d Cohen's d
#' @param n1 Sample size in group 1
#' @param n2 Sample size in group 2
#' 
#' @details See page 7 in Erdfelder, Faul, & Buchner (1996) for the
#'     formula.
#' 
#' @return The non-centrality parameter
#' 
#' @references 
#' Erdfelder, E., Faul, F., & Buchner, A. (1996). GPOWER: A general 
#' power analysis program. Behavior research methods, instruments, & 
#' computers, 28(1), 1-11.
#'
#' @examples
#' ncp_from_d(0.3, 50, 50)
#'
#' @export
#' 
ncp_from_d <- function(d, n1, n2) {
  ncp <- d * sqrt((n1 * n2^2 + n2 * n1^2) / ((n1 + n2)^2))
  return(ncp)
}
