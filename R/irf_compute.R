#' Compute global structural IRFs (VAR(1))
#'
#' @description
#' Computes global-level structural impulse responses from posterior
#' draws of the global VAR parameters (\code{A} and \code{Sigma}).
#'
#' @param chains A list of MCMC chains from \code{run_mcmc()$CHAINS}.
#' @param H Number of horizons.
#'
#' @return A 4D array:
#' \code{Phi_global[draw, horizon, response, shock]}.
#'
#' @export
compute_global_irf <- function(chains, H){

    cc <- combine_global_chains(chains)

    Nsim_tot <- dim(cc$A)[1]
    N        <- nrow(cc$S[1,,])

    Phi_global <- array(NA_real_, dim = c(Nsim_tot, H+1, N, N))

    for (s in 1:Nsim_tot){

        A <- cc$A[s,,]
        S <- cc$S[s,,]
        S <- (S + t(S))/2
        P <- chol(S + diag(1e-8, N))

        psi <- diag(N)
        Phi_global[s,1,,] <- psi %*% P

        for (h in 2:(H+1)){
            psi <- psi %*% A
            Phi_global[s,h,,] <- psi %*% P
        }
    }

    return(Phi_global)
}

