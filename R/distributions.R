# distributions.R

# Generators with additional numerical stability guaranteed by:
# 1. Systematic symmetrization
# 2. Ridge coefficient
# 3. Cholesky-only inversion

#' Matrix-normal generator
#'
#' @description
#' Generates a draw from a matrix-normal distribution
#' with mean matrix \code{M}, row covariance \code{V}, and
#' column covariance \code{S}.
#'
#' @param M Matrix. Mean matrix of size K × N.
#' @param V Matrix. Row covariance matrix (K × K), positive-definite.
#' @param S Matrix. Column covariance matrix (N × N), positive-definite.
#'
#' @return A numeric matrix of dimension K × N.
#'
#' @examples
#' M  <- matrix(0, 3, 2)
#' V  <- diag(3)
#' S  <- diag(2)
#' rmatnorm(M, V, S)
#'
#' @export
rmatnorm <- function(M, V, S) {
    M <- as.matrix(M); V <- as.matrix(V); S <- as.matrix(S)
    K <- nrow(M); N <- ncol(M)

    Uv <- chol(V + diag(1e-8, nrow(V)))
    Us <- chol(S + diag(1e-8, nrow(S)))

    Z <- matrix(rnorm(K * N), K, N)
    A <- M + t(Uv) %*% Z %*% Us
    return(A)
}

#' Wishart generator
#'
#' @description
#' Generates a random draw from a Wishart distribution with
#' degrees of freedom \code{df} and scale matrix \code{S}.
#'
#' @param df Integer. Degrees of freedom (df ≥ dimension of S).
#' @param S Matrix. Scale matrix (K × K), symmetric and positive-definite.
#'
#' @return A numeric positive-definite matrix drawn from W(df, S).
#'
#' @examples
#' S <- diag(3)
#' wishart_rnd(10, S)
#'
#' @export
wishart_rnd <- function(df, S) {
    S <- (S + t(S))/2
    K <- nrow(S)
    S <- S + diag(1e-8, K)
    W <- stats::rWishart(1, df, S)[,,1]
    W <- (W + t(W))/2
    W <- W +diag(1e-8,nrow(W))

}

#' Inverse-Wishart generator
#'
#' @description
#' Generates a draw from an inverse-Wishart distribution with
#' degrees of freedom \code{df} and scale matrix \code{S}.
#'
#' @param df Integer. Degrees of freedom (df ≥ dimension of S).
#' @param S Matrix. Scale matrix (K × K), symmetric and positive-definite.
#'
#' @return A numeric positive-definite matrix drawn from IW(df, S).
#'
#' @examples
#' S <- diag(3)
#' invwishart_rnd(10, S)
#'
#' @export
invwishart_rnd <- function(df, S) {
    S <- (S + t(S))/2
    K <- nrow(S)
    S <- S + diag(1e-8, K)
    Winv_scale <- chol2inv(chol(S))
    W <- stats::rWishart(1, df, Winv_scale)[,,1]
    U <- chol(W + diag(1e-8,nrow(W)))
    Wi <- chol2inv(U)
    (Wi + t(Wi))/2
}
