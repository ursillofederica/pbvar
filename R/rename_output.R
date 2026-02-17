#' Add names to MCMC output (internal)
#'
#' @description
#' Internal helper used by \code{run_mcmc()} to assign row and column
#' names to all sampled matrices in the MCMC output.
#' It attaches:
#' \itemize{
#'   \item variable names to rows of \code{A}, \code{Ag}, \code{V}, \code{Sigma}, \code{Sg};
#'   \item lag names to columns of \code{A} and \code{Ag};
#'   \item country names to the column dimension of \code{Ag} and \code{Sg};
#'   \item reshapes \code{Ag} and \code{Sg} into matrices of size
#'         \code{eff_iter × n_groups}.
#' }
#'
#'
#' @param chains A list of MCMC chains produced by \code{run_mcmc()}.
#' @param var_names Character vector of variable names.
#' @param lag_names Character vector of lag names.
#' @param country_names Character vector of country names.
#' @param n_groups Integer, number of countries.
#' @param R Total number of post–burn-in iterations (before thinning).
#' @param thin Thinning interval used in the sampler.
#'
#' @return The same \code{chains} object, but with added dimension names.
#'
#' @keywords internal
#' @noRd
rename_output <- function(chains, var_names, lag_names, country_names, n_groups, R, thin){

    eff_iter <- as.integer(R/thin)

    for (ch in 1:length(chains)) {

        for (i in seq_along(chains[[ch]]$samples$A)) {
            M <- chains[[ch]]$samples$A[[i]]
            rownames(M) <- var_names
            colnames(M) <- lag_names
            chains[[ch]]$samples$A[[i]] <- M
        }

        for (i in seq_along(chains[[ch]]$samples$V)) {
            M <- chains[[ch]]$samples$V[[i]]
            rownames(M) <- var_names
            colnames(M) <- var_names
            chains[[ch]]$samples$V[[i]] <- M
        }

        for (i in seq_along(chains[[ch]]$samples$Sigma)) {
            M <- chains[[ch]]$samples$Sigma[[i]]
            rownames(M) <- var_names
            colnames(M) <- var_names
            chains[[ch]]$samples$Sigma[[i]] <- M
        }
    }

    for (ch in 1:length(chains)) {
        dim(chains[[ch]]$samples$Ag) <- c(eff_iter, n_groups)
        dim(chains[[ch]]$samples$Sg) <- c(eff_iter, n_groups)

        colnames(chains[[ch]]$samples$Ag) <- country_names
        colnames(chains[[ch]]$samples$Sg) <- country_names
    }

    for (ch in 1:length(chains)) {
        for (j in 1:n_groups) {
            for (i in 1:eff_iter) {

                M <- chains[[ch]]$samples$Ag[[i, j]]
                rownames(M) <- var_names
                colnames(M) <- lag_names
                chains[[ch]]$samples$Ag[[i, j]] <- M

                S <- chains[[ch]]$samples$Sg[[i, j]]
                rownames(S) <- var_names
                colnames(S) <- var_names
                chains[[ch]]$samples$Sg[[i, j]] <- S
            }
        }
    }
    return(chains)
}
