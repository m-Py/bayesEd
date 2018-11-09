
#' Compute "power" for the default t-test Bayes factor
#'
#' For a given N, this function returns the Bayes factor that will be
#' exceeded with a specified probability in a t-test. Transfers the
#' frequentist concept of statistical power to Bayes factors. Uses the
#' default Bayesian t test implemented in the package `BayesFactor`.
#' 
#' @param N The total sample size over both groups in the independent
#'     groups t test.
#' @param effect_size The assumed effect size d.
#' @param nsim The number of simulations that are conducted to determine
#'     the "power" of the t test. Increase this number to obtain a more
#'     precise estimate.
#' @param rscale The scaling parameter in the Cauchy prior used in the
#'     Bayes factor computation. Defaults to `sqrt(2) / 2`.
#' @param probability A scalar indicating the "power". See `Details`.
#' @param say_result Boolean. Print a verbal description of the
#'     computation when it is finished?
#' 
#' @details For Bayes factors, there is no concept of statistical
#'     power. However, for a given N and effect size, we can investigate
#'     the distribution of Bayes factors. Here, "power" then means the
#'     probability that an observed Bayes factor is at least as high as
#'     a particular value. We specify this "power" using the parameter
#'     `probability`. In case of a null-effect, "power" is the
#'     probability that an observed Bayes factor is smaller than a
#'     particular value. The function does not analytically determine
#'     "power", but uses simulation to generate the distribution of
#'     Bayes factors.
#' 
#' @examples
#' ## N is the total sample size across two groups
#' power_bf(N = 200, effect_size = 0.4)
#' ## Null effect:
#' power_bf(N = 1000, effect_size = 0)
#' ## Small effect size:
#' power_bf(N = 1000, effect_size = 0.2)
#' ## Vary width of the prior distribution for the alternative hypothesis:
#' power_bf(N = 1000, effect_size = 0.2, rscale = 0.2)
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

power_bf <- function(N, effect_size, nsim = 1000, rscale = sqrt(2)/2,
                     probability = 0.8, say_result = TRUE) {
  if (effect_size > 0) {
    power <- 1 - probability
    direction <- "larger"
  } else if (effect_size == 0) {
    power <- probability
    direction <- "smaller"
  } else {
    stop("effect size cannot be lower than 0")
  }
  info <- BF_quantiles(estimate_expected_bf(effect_size, N, nsim, rscale), power)
  info$BF <- exp(info$logBF) # other functions return the log BF
  info$quantile <- NULL
  info$power    <- probability
  if (say_result) {
    cat("\n Approximately", probability*100,
        "% of all Bayes factors will be", direction, "than",
        round(info$BF, 2), "\n\n")
  }
  return(info)
}

#' Simulate the distribution of the t-test Bayes factor
#'
#' Uses the default Bayes factor employed in the BayesFactor package.
#'
#' @param effect_size The assumed effect size d.
#' @param sample_sizes A vector of sample sizes to be used in the
#'     simulation. Means the total sample sizes, i.e., the sample sizes
#'     across the two groups in the t test.
#' @param n_bayes_factors How many Bayes factor should be computed. Is
#'     replicated for each of the elements in the vector `sample_sizes`.
#' @param rscale The scaling parameter in the Cauchy prior used in the
#'     Bayes factor computation. Defaults to `sqrt(2) / 2`.
#'
#' @return A data.frame in long format where each row represents the
#'     simulation of a t test. Each row has four columns: `logBF` - the
#'     natural logarithm of the Bayes factor; `N` - the total sample
#'     size in the t test; `eff_size` - the effect size; `rscale` - the
#'     scaling parameter of the prior distribution.
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
#' @examples 
#' bfs <- bf_distribution(effect_size = 0.5, sample_sizes = seq(50, 350, by = 50), n_bayes_factors = 300)
#'
#' ## Use `exp` to obtain the "normal" BF (not logarithm of the Bayes factor)
#' tapply(exp(bfs$logBF), bfs$N, median)
#'
#' 
#' @export
#' 
#' @importFrom BayesFactor ttestBF
#' @importFrom BayesFactor extractBF
#' 
bf_distribution <- function(effect_size, sample_sizes,
                            n_bayes_factors, rscale=sqrt(2)/2) {
  ## store BFs for each repetition for one sample size:
  repetitions_bf <- vector(length = n_bayes_factors)
  ## store BFs by sample size:
  samples_bf     <- list()
  for (i in 1:length(sample_sizes)) {
    for (j in 1:n_bayes_factors) {
      group0 <- rnorm(sample_sizes[i]/2, 0, 1)
      group1 <- rnorm(sample_sizes[i]/2, effect_size, 1)
      tmpBF  <- BayesFactor::ttestBF(group0, group1, rscale = rscale)
      repetitions_bf[j] <- BayesFactor::extractBF(tmpBF)$bf
    }
    ## store all bayes factors for sample size i (store the log BF!)
    samples_bf[[i]] <- log(repetitions_bf)
  }
  ## convert to data.frame in long format:
  ret <- data.frame(logBF = unlist(samples_bf),
                    N  = rep(sample_sizes, each = n_bayes_factors),
                    eff_size = effect_size,
                    rscale = rscale)
  return(ret)
}

#' Compute quantiles for the distribution of Bayes factors
#'
#' @param dat A `data.frame` returned by `estimate_expected_bf`
#' @param quantiales The quantiles computed for the distribution of
#'     Bayes factors. Passed to function `quantile`.
#' @return A `data.frame` in long format containing Bayes factor
#'     quantiles per sample size.
#' @importFrom reshape2 melt
#' @importFrom plyr adply
#'
#' @examples
#' 
#' bfs <- bf_distribution(effect_size = 0.5, sample_sizes = seq(50, 350, by = 50), n_bayes_factors = 300)
#' BF_quantiles(bfs)
#'
#' @export
#' 
BF_quantiles <- function(dat, quantiles = c(0.2, 0.5, 0.8)) {
  ## obtain quantiles of log BF by sample size
  arr <- tapply(dat$logBF, dat$N, quantile, probs = quantiles,
                simplify = FALSE)
  ## use plyr to convert array to data.frame
  dfw <- plyr::adply(arr, 1)
  ## use reshape to create long data
  long_quantiles <- reshape2::melt(dfw, id.vars = c("X1"), factorsAsStrings = TRUE)
  
  ## Some renaming:
  long_quantiles$N <- as.numeric(levels(long_quantiles$X1))[long_quantiles$X1]
  long_quantiles$X1 <- NULL
  long_quantiles$quantile <- long_quantiles$variable
  long_quantiles$variable <- NULL
  long_quantiles$logBF    <- long_quantiles$value
  long_quantiles$value    <- NULL
  
  ## Attach effect size and prior scale
  long_quantiles$eff_size <- unique(dat$eff_size)
  long_quantiles$rscale   <- unique(dat$rscale)
  return(long_quantiles)
}

#' Plot the distribution of Bayes factors
#'
#' Plots  quantiles of Bayes factors by sample size.
#'
#' @param BF_quantiles A `data.frame` returned by `BF_quantiles`.
#' @param thresholds A vector of thresholds to be drawn
#'     horizontally. May for example illustrate some thresholds over
#'     which Bayes factors are deemed to be informative.
#' @param ylim A vector of two values indicating the limits of the
#'     y-axis. Given on a log-scale.
#' @param main The title of the plot.
#' @param axis A character vector indicating what y-axis is to
#'     drawn. Can be "log" for the log-axis, "standard" for
#'     non-logarithmic axis, or both (i.e., c("log", "standard") which
#'     is also the default)
#' 
#' @examples 
#' bfs <- bf_distribution(effect_size = 0.5, sample_sizes = seq(50, 350, by = 50), n_bayes_factors = 300)
#' quants <- BF_quantiles(bfs)
#' plot_BF_quantiles(quants, thresholds = c()
#' legend("topleft", legend = paste0("PR ", c(80, 50, 20)), pch = 3, col = 3:1, bty = "n")
#'
#' @export

plot_BF_quantiles <- function(BF_quantiles, thresholds=NULL, ylim=NULL,
                              main = "", axis = c("log", "standard")) {
  
  ## first make some adjustments to the plot in dependence of which
  ## axis the user wants to show:
  ylab <- ""
  ## adjust margin to make space for a second y-axis -- ‘c(bottom,
  ## left, top, right)’; The default ‘c(5, 4, 4, 2) + 0.1’.
  def_mar <- c(5, 4, 4, 2) + 0.1
  cur_mar <- def_mar
  if ("standard" %in% axis) {
    cur_mar <- cur_mar + c(0,0,0,2)
    par(mar = cur_mar)
  }
  else {
    cur_mar <- cur_mar + c(0,0,0,-2)
    par(mar = cur_mar)
  }
  if ("log" %in% axis) {
    draw_y <- "s"
    ylab <- "log(BF10)"
  } else {
    draw_y <- "n"
    cur_mar <- cur_mar + c(0,-2,0,0)
    par(mar = cur_mar)
  }
  
  ## draw the actual plot
  plot(BF_quantiles$N, BF_quantiles$logBF, pch = 3,
       col = rep(1:length(unique(BF_quantiles$quantile)),
                 each = length(unique(BF_quantiles$N))), ylab = ylab,
       xlab = "", las = 1, ylim = ylim, main = main,
       yaxt = draw_y)
  
  ## add y-axis to the right if that was required by the user
  if ("standard" %in% axis) 
    add_normal_bf_scale(1)
  
  ## indicate the neutral evidence line if it was required by the user
  cols = rep("darkgrey", length(thresholds))
  cols[which(thresholds == 1)] <- "black"
  if (!is.null(thresholds))
    abline(h = log(thresholds), lty = 2, lwd = 1.5, col = cols)
  
  ## Label x-axis
  axis(1, 20, "n", tick = FALSE)
  
  ## reset margins
  on.exit(par(mar = def_mar))
}

## Add an axis on the right side of the plot for "normal" BFs in
## addition to log BFs
add_normal_bf_scale <- function(by = 1) {
  label_bf_lim <- seq(-30, 30, by = by)
  labs <- round(exp(label_bf_lim), 2)
  labs[labs >= 10] <- round(labs[labs >= 10])
  labs[labs >= 1] <- round(labs[labs >= 1], 1)
  axis(4, at = label_bf_lim, labels = labs, las = 1)
  mtext("BF10", 4, line = 3.1)
}
