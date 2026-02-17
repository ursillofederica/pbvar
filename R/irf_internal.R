#' Combine Ag and Sg draws across chains
#'
#' @description
#' Internal helper that vertically stacks country-specific coefficient
#' and covariance draws from multiple MCMC chains. Each chain stores
#' \code{Ag} and \code{Sg} as matrices of lists, with dimension
#' \code{(iterations × G)}. This function binds them across chains.
#'
#' @param chains A list of chain outputs as stored inside
#'   \code{run_mcmc()$CHAINS}.
#'
#' @return
#' A list with:
#' \itemize{
#'   \item \code{Ag}: matrix of combined coefficient draws (K × N),
#'   \item \code{Sg}: matrix of combined covariance draws (N × N).
#' }
#'
#' @keywords internal
#' @noRd
combine_group_chains <- function(chains){

    nChains <- length(chains)
    Ag_all <- chains[[1]]$samples$Ag
    Sg_all <- chains[[1]]$samples$Sg

    for( i in 2:nChains){

        Ag_all <- rbind(Ag_all, chains[[i]]$samples$Ag)
        Sg_all <- rbind(Sg_all, chains[[i]]$samples$Sg)

    }

    list(Ag = Ag_all, Sg = Sg_all)
}

#' Posterior IRF summary statistics
#'
#' @description
#' Computes posterior median and 95% equal-tailed intervals for the
#' structural impulse responses obtained from
#' \code{compute_group_irf()}.
#'
#' For each country, impulse variable, response variable, and horizon,
#' the function extracts:
#' \itemize{
#'   \item posterior median,
#'   \item 2.5% quantile,
#'   \item 97.5% quantile.
#' }
#'
#' @param Phi A 5D array returned by \code{compute_group_irf()}.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{phi_med}: posterior medians,
#'   \item \code{phi_low}: lower bound,
#'   \item \code{phi_high}: upper bound.
#' }
#'
#' @keywords internal
#' @noRd
stats_group_irf <- function(Phi){

    Nsim <- dim(Phi)[1]
    G    <- dim(Phi)[2]
    H1   <- dim(Phi)[3]
    N    <- dim(Phi)[4]

    phi_med  <- array(NA, dim = c(G, H1, N, N))
    phi_low  <- array(NA, dim = c(G, H1, N, N))
    phi_high <- array(NA, dim = c(G, H1, N, N))

    for (g in 1:G){
        for (h in 1:H1){
            for (i in 1:N){
                for (j in 1:N){


                    draws_ij            <- Phi[, g, h, i, j]
                    phi_med[g, h, i, j] <- median(draws_ij)

                    q                    <- quantile(draws_ij, probs = c(0.025, 0.975))
                    phi_low[g, h, i, j]  <- q[1]
                    phi_high[g, h, i, j] <- q[2]
                }
            }
        }
    }

    return(list(phi_med  = phi_med, phi_low  = phi_low, phi_high = phi_high))
}

#' Combine global-level draws across MCMC chains (internal)
#'
#' @description
#' Internal helper that stacks posterior draws of the global VAR
#' parameters \code{A} and \code{Sigma} across multiple MCMC chains.
#' Each chain stores \code{A} and \code{Sigma} as lists of matrices.
#'
#' @param chains A list of chains returned by \code{run_mcmc()$CHAINS}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{A}: 3D array \code{[draw × K × N]} of coefficient matrices,
#'   \item \code{S}: 3D array \code{[draw × N × N]} of covariance matrices.
#' }
#'
#' @keywords internal
#' @noRd
combine_global_chains <- function(chains){

    nChains <- length(chains)
    nEff    <- length(chains[[1]]$samples$A)

    A0     <- chains[[1]]$samples$A[[1]]
    Sigma0 <- chains[[1]]$samples$Sigma[[1]]

    K <- nrow(A0)
    N <- ncol(A0)
    Nsim_tot <- nChains * nEff

    A_all     <- array(NA_real_, dim = c(Nsim_tot, K, N))
    Sigma_all <- array(NA_real_, dim = c(Nsim_tot, N, N))

    idx <- 0
    for (c in 1:nChains) {
        for (r in 1:nEff) {
            idx <- idx + 1
            A_all[idx,,]     <- chains[[c]]$samples$A[[r]]
            Sigma_all[idx,,] <- chains[[c]]$samples$Sigma[[r]]
        }
    }

    return(list(A=A_all, S=Sigma_all))
}
