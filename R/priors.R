#' Sample hyperparameters for the VAR model
#'
#' @description
#' Draws a set of scalar hyperparameters from their respective
#' prior distributions, using the values stored in the list
#' \code{STV}.
#' Each element is generated from a univariate distribution:
#' normal or gamma/inverse-gamma, depending on the parameter.
#'
#' @param STV List containing the hyperparameter values used as
#' inputs for the prior distributions. Expected elements include
#' \code{mum}, \code{sigma2m}, \code{sw}, \code{aw}, \code{ss},
#' \code{nus}, \code{as}, and \code{bs}.
#'
#' @return
#' A list with sampled hyperparameters:
#' \itemize{
#'   \item \code{m} — normal draw with mean \code{mum} and variance \code{sigma2m}.
#'   \item \code{s} — inverse-gamma-like draw: \eqn{1 / \mathrm{Gamma}(ss,\ nus)}.
#'   \item \code{w} — gamma draw with shape \code{sw} and scale \code{aw}.
#'   \item \code{s_Sigma} — inverse-gamma-like draw: \eqn{1 / \mathrm{Gamma}(as,\ bs)}.
#' }
#'
#' @examples
#' STV <- build_starting_values(N = 3, K = 6, lag = 2)
#' hyp <- hyper(STV)
#' str(hyp)
#'
#' @export
hyper <- function(STV) {
    mum     <- STV$mum
    sigma2m <- STV$sigma2m
    sw      <- STV$sw
    aw      <- STV$aw
    ss      <- STV$ss
    nus     <- STV$nus
    as      <- STV$as
    bs      <- STV$bs

    H <- list()
    H$m       <- rnorm(1, mean = mum, sd = sqrt(sigma2m))
    H$s       <- 1 / rgamma(1, shape = ss, scale = nus)
    H$w       <- rgamma(1, shape = sw, scale = aw)
    H$s_Sigma <- 1 / rgamma(1, shape = as, scale = bs)

    H
}
