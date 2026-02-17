#' Posterior summaries for panel VAR parameters
#'
#' @description
#' Computes posterior summary statistics (mean, sd, quantiles, HPD interval,
#' Gelman–Rubin R-hat, and effective sample size) for selected parameter blocks
#' of a hierarchical panel VAR model estimated with \code{run_mcmc()}.
#'
#' The function can summarize:
#' \itemize{
#'   \item \code{A}     — global autoregressive coefficients;
#'   \item \code{Ag}    — country-specific autoregressive coefficients;
#'   \item \code{Sg}    — country-specific covariance matrices;
#'   \item \code{V}     — global coefficient covariance matrix;
#'   \item \code{Sigma} — global error covariance matrix.
#' }
#'
#' The user can select which blocks to summarize via the \code{what} argument.
#'
#' @param chains A list of chains as returned in \code{run_mcmc()$CHAINS}.
#' @param var_names Character vector with variable names (row names of A and Sg).
#' @param lag_names Character vector with lag names (column names of A and Ag).
#' @param country_names Character vector with country identifiers (columns of Ag and Sg).
#' @param par Character vector indicating which parameters to summarize.
#'   Possible values: \code{"A"}, \code{"Ag"}, \code{"Sg"},
#'   \code{"V"}, \code{"Sigma"}, or \code{"all"}.
#'
#' @return A named list of data frames, one for each selected parameter block.
#'
#'
#' @examples
#' # Minimal example using a single chain
#' Y_list <- list(matrix(rnorm(50), 10, 5))
#' Z_list <- list(matrix(rnorm(50), 10, 5))
#'
#' out <- run_mcmc(
#'   Y_list, Z_list,
#'   lag = 1,
#'   fixed_list = NULL,
#'   R = 100,
#'   burnin = 50,
#'   thin = 2,
#'   nChains = 2,
#'   var_names = paste0("V", 1:5),
#'   lag_names = paste0("L", 1:5),
#'   country_names = "C1"
#' )
#'
#' summary_par(out$CHAINS,
#'             var_names = paste0("V", 1:5),
#'             lag_names = paste0("L", 1:5),
#'             country_names = "C1",
#'             par = "A")
#' @export
summary_par <- function(chains,var_names,lag_names,country_names, par = c("A","Ag","Sg","V","Sigma","all")){

    summary_A <- list()

    for (i in 1:3) {
        for (j in 1:3) {

            mlist <- coda::mcmc.list(
                lapply(seq_along(chains), function(cc)
                    coda::mcmc(sapply(chains[[cc]]$samples$A, function(M) M[i,j]))
                )
            )
            name <- paste0("A[", var_names[i], ",", lag_names[j], "]")
            summary_A[[name]] <- summarize_mcmc(mlist)
        }
    }

    summary_Ag <- list()

    for (country in country_names) {
        for (i in 1:3) {
            for (j in 1:3) {

                mlist <- coda::mcmc.list(
                    lapply(seq_along(chains), function(cc)
                        coda::mcmc(sapply(chains[[cc]]$samples$Ag[, country], function(M) M[i,j]))
                    )
                )

                name <- paste0("Ag[", country, ", ", var_names[i], ", ", lag_names[j], "]")
                summary_Ag[[name]] <- summarize_mcmc(mlist)
            }
        }
    }

    summary_Sg <- list()

    for (country in country_names) {
        for (i in 1:3) {
            for (j in 1:3) {
                mlist <- coda::mcmc.list(
                    lapply(seq_along(chains), function(cc)
                        coda::mcmc(sapply(chains[[cc]]$samples$Sg[, country], function(M) M[i,j]))
                    )
                )
                name <- paste0("Sg[", country, ", ", var_names[i], ", ", var_names[j], "]")
                summary_Sg[[name]] <- summarize_mcmc(mlist)
            }
        }
    }



    summary_V <- list()

    for (i in 1:3) {
        for (j in 1:3) {

            mlist <- coda::mcmc.list(
                lapply(seq_along(chains), function(cc)
                    coda::mcmc(sapply(chains[[cc]]$samples$V, function(M) M[i,j]))
                )
            )

            name <- paste0("V[", var_names[i], ", ", var_names[j], "]")
            summary_V[[name]] <- summarize_mcmc(mlist)
        }
    }

    summary_Sigma <- list()

    for (i in 1:3) {
        for (j in 1:3) {

            mlist <- coda::mcmc.list(
                lapply(seq_along(chains), function(cc)
                    coda::mcmc(sapply(chains[[cc]]$samples$Sigma, function(M) M[i,j]))
                )
            )
            name <- paste0("Sigma[", var_names[i], ", ", var_names[j], "]")
            summary_Sigma[[name]] <- summarize_mcmc(mlist)
        }
    }

    tidy_A     <- tidy_summary(summary_A)
    tidy_V     <- tidy_summary(summary_V)
    tidy_Sigma <- tidy_summary(summary_Sigma)
    tidy_Ag    <- tidy_summary(summary_Ag)
    tidy_Sg    <- tidy_summary(summary_Sg)

    results <- list()

    if ("A"     %in% par || "all" %in% par) results$A     <- tidy_A
    if ("Ag"    %in% par || "all" %in% par) results$Ag    <- tidy_Ag
    if ("Sg"    %in% par || "all" %in% par) results$Sg    <- tidy_Sg
    if ("V"     %in% par || "all" %in% par) results$V     <- tidy_V
    if ("Sigma" %in% par || "all" %in% par) results$Sigma <- tidy_Sigma

    return(results)
}
