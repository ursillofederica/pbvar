#' Build starting values for the sampler
#'
#' @description
#' Creates a list of initial parameter values and hyperparameters
#' needed to start the sampler. The function builds
#' defaults for matrices \code{M}, \code{S}, \code{SS}, \code{W},
#' and scalar hyperparameters, optionally replacing some
#' of them with user-provided fixed values.
#'
#' @param N Integer. Number of  variables.
#' @param K Integer. State dimension of the companion form (K = N * lag)
#' of the state matrix.
#' @param lag Integer. Maximum autoregressive lag of the VAR.
#' @param fixed_list Optional list. If provided, it can contain any
#' of the elements \code{m}, \code{sigma2m}, \code{sw}, \code{aw},
#' \code{ss}, \code{nus}, \code{as}, \code{bs}. Each non-NULL entry
#' replaces the corresponding default value in the output.
#'
#' @return
#' A list containing starting values for:
#' \itemize{
#'   \item \code{M}  – K × N mean matrix for global parameter A.
#'   \item \code{S}, \code{SS} – N × N covariance-related matrices.
#'   \item \code{W} – K × K diagonal matrix based on decreasing weights based on lags.
#'   \item fixed degrees of freedom: \code{mu_Sigma}, \code{eta}, \code{nu} chosen to ensure finite
#'   first moment for inverse-Wishart and Wishart distributions.
#'   And the elements possibly taken from \code{fixed_list}.
#' }
#'
#' @examples
#' init <- build_starting_values(N = 5, K = 10, lag = 2)
#' @export
build_starting_values <- function(N, K, lag, fixed_list = NULL) {
    STV <- list()

    STV$M  <- matrix(0, nrow = K, ncol = N)
    v      <- numeric(K)
    STV$S  <- diag(1, N)
    STV$SS <- diag(1, N)

    row <- 1L
    for (i in seq_len(lag)) {
        for (j in seq_len(N)) {
            if (i == 1L) {
                STV$M[row, j] <- 1
            }
            v[row] <- 1 / (i^2)
            row <- row + 1L
            if (row > K) break
        }
        if (row > K) break
    }

    STV$W <- diag(v, nrow = K, ncol = K)

    STV$mu_Sigma <- N + 10
    STV$eta      <- K + 8
    STV$nu       <- N + 6

    STV$mum      <- if (!is.null(fixed_list$m))        fixed_list$m        else 0
    STV$sigma2m  <- if (!is.null(fixed_list$sigma2m))  fixed_list$sigma2m  else 1
    STV$sw       <- if (!is.null(fixed_list$sw))       fixed_list$sw       else 2
    STV$aw       <- if (!is.null(fixed_list$aw))       fixed_list$aw       else 1
    STV$ss       <- if (!is.null(fixed_list$ss))       fixed_list$ss       else 2
    STV$nus      <- if (!is.null(fixed_list$nus))      fixed_list$nus      else 0.5
    STV$as       <- if (!is.null(fixed_list$as))       fixed_list$as       else 2
    STV$bs       <- if (!is.null(fixed_list$bs))       fixed_list$bs       else 0.5


    STV
}

#' Build initial state matrices for the sampler
#'
#' @description
#' Computes starting values for the state of the sampler in a panel-VAR
#' or multi-group VAR setting. For each group, the function estimates
#' initial coefficient matrices using OLS and computes corresponding
#' residual covariance matrices. It then combines the group-specific
#' quantities to produce initial values for \code{A}. \code{Sigma_g}, \code{Sigma},
#' and \code{V} based on MAP estimates under Wishart and inverse-Wishart priors.
#'
#' @param STV List of starting values produced by
#' \code{build_starting_values()}, containing \code{M}, \code{W},
#' \code{SS}, \code{mu_Sigma}, \code{eta}, \code{nu}, and related
#' hyperparameters.
#' @param Y_list List of data matrices \code{Y_i}, each of size
#' \eqn{T_i \times N}, one for each group.
#' @param Z_list List of regressor matrices \code{Z_i}, each of size
#' \eqn{T_i \times K}, one for each group, typically the stacked-lag
#' companion regressors.
#' @param Hyper List of scalar hyperparameters sampled by \code{hyper()},
#' including \code{m}, \code{s}, \code{w}, and \code{s_Sigma}.
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item \code{A} – K × N pooled initial coefficient matrix.
#'   \item \code{Sigma} – N × N pooled residual covariance matrix
#'   (MAP inverse-Wishart estimate).
#'   \item \code{V} – K × K state covariance matrix (MAP inverse-Wishart estimate).
#'   \item \code{Ag} – list of group-specific OLS coefficient matrices (K × N).
#'   \item \code{Sg} – list of group-specific covariance matrices (N × N).
#' }
#'
#' @details
#' For the covariance matrices, the function uses MAP estimates under
#' inverse-Wishart/Wishart priors instead of direct prior draws. This
#' provides numerically stable and reasonably central starting values
#' for the sampler, which is particularly useful in high-dimensional
#' settings with many parameters.
#'
#' \enumerate{
#'   \item Group-wise OLS estimation of VAR coefficients using a QR decomposition.
#'   \item Computation of group residual covariance matrices \code{RSS_i}.
#'   \item Pooled MAP estimates for \code{Sigma} and \code{V} under inverse-Wishart priors.
#'   \item Group-specific MAP estimates \code{Sg} for residual covariances.
#' }
#'
#' @examples
#' # Example with dummy data (small dimensions)
#' Y_list <- list(matrix(rnorm(50), 10, 5))
#' Z_list <- list(matrix(rnorm(100), 10, 10))
#' STV <- build_starting_values(N = 5, K = 10, lag = 2)
#' Hyper <- hyper(STV)
#' out <- init(STV, Y_list, Z_list, Hyper)
#' str(out)
#'
#' @export
init <- function(STV, Y_list, Z_list, Hyper) {
    G <- length(Y_list)
    N <- ncol(Y_list[[1]])
    K <- ncol(Z_list[[1]])

    Ag_init <- vector("list", G)
    Sg_init <- vector("list", G)
    RSSg    <- vector("list", G)

    T_tot   <- 0L
    RSS_tot <- matrix(0, N, N)

    # OLS per Ag_init
    for (i in seq_len(G)) {
        Zi <- Z_list[[i]]
        Yi <- Y_list[[i]]

        qrZ <- qr(Zi)
        Q   <- qr.Q(qrZ)
        R   <- qr.R(qrZ)

        QtY <- crossprod(Q, Yi)
        Agi <- backsolve(R, QtY)

        Ag_init[[i]] <- Agi              # K x N

        Ei      <- Yi - Zi %*% Agi       # T_i x N
        RSSi    <- crossprod(Ei)         # N x N
        RSSi    <- (RSSi + t(RSSi)) / 2
        RSSg[[i]] <- RSSi

        RSS_tot <- RSS_tot + RSSi
        T_tot   <- T_tot + nrow(Yi)
    }

    # Media iniziale di A come media di Ag_init
    A_stack <- array(NA_real_, dim = c(K, N, G))
    for (i in seq_len(G)) A_stack[,,i] <- Ag_init[[i]]
    A_init <- apply(A_stack, c(1, 2), mean)   # K x N

    # Sigma iniziale (MAP IW)
    RSS_tot       <- (RSS_tot + t(RSS_tot)) / 2
    S_Sigma_post  <- Hyper$s * STV$SS + RSS_tot      # N x N
    S_Sigma_post  <- (S_Sigma_post + t(S_Sigma_post)) / 2
    df_S_post     <- STV$mu_Sigma + T_tot

    Sigma_init <- S_Sigma_post / (df_S_post + N + 1) # MAP IW
    Sigma_init <- (Sigma_init + t(Sigma_init)) / 2

    # Sg (MAP IW)
    RSSa_tot <- matrix(0, K, K)
    for (i in seq_len(G)) {
        Ea     <- Ag_init[[i]] - A_init           # K x N
        RSSa   <- Ea %*% t(Ea)                    # K x K
        RSSa   <- (RSSa + t(RSSa)) / 2
        RSSa_tot <- RSSa_tot + RSSa

        S_Sigmag_post  <- (STV$nu - N - 1) * Sigma_init + RSSg[[i]]  # N x N
        S_Sigmag_post  <- (S_Sigmag_post + t(S_Sigmag_post)) / 2
        nu_Sigmag_post <- STV$nu + nrow(Y_list[[i]])

        Sg_i <- S_Sigmag_post / (nu_Sigmag_post + N + 1)  # MAP IW
        Sg_init[[i]] <- (Sg_i + t(Sg_i)) / 2
    }

    # V iniziale (MAP IW)
    S_V_post  <- Hyper$w * STV$W + RSSa_tot     # K x K
    S_V_post  <- (S_V_post + t(S_V_post)) / 2
    df_V_post <- STV$eta + G * N

    V_init <- S_V_post / (df_V_post + K + 1)
    V_init <- (V_init + t(V_init)) / 2
    V_init <- V_init + diag(1e-8, nrow(V_init))

    State <- list()
    State$A     <- A_init       # K x N
    State$V     <- V_init       # K x K
    State$Sigma <- Sigma_init   # N x N

    State$Ag <- vector("list", G)
    State$Sg <- vector("list", G)
    for (g in seq_len(G)) {
        State$Ag[[g]] <- Ag_init[[g]]  # K x N
        State$Sg[[g]] <- Sg_init[[g]]  # N x N
    }

    State
}
