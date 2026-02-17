#' Run the hierarchical VAR MCMC sampler
#'
#' @description
#' Runs the full MCMC algorithm for the hierarchical VAR model.
#' For each chain, the function:
#' \enumerate{
#'   \item builds starting values and initializes hyperparameters;
#'   \item performs the group-specific MNIW updates;
#'   \item updates global parameters \code{A}, \code{V}, and \code{Sigma};
#'   \item updates all scalar hyperparameters;
#'   \item stores thinned post-burnin draws.
#' }
#' Multiple chains can be run independently, each with its own seed.
#'
#' @param Y_list List of matrices \code{Y_g}, each of size \eqn{T_g \times N},
#' representing the dependent variables for each group.
#' @param Z_list List of matrices \code{Z_g}, each of size \eqn{T_g \times K},
#' containing the stacked-lag regressors of the VAR.
#' @param lag Integer. Maximum autoregressive lag of the VAR.
#' @param fixed_list Optional list of fixed hyperparameters passed to
#' \code{build_starting_values()}.
#' @param R Integer. Number of post-burnin MCMC iterations (before thinning).
#' @param burnin Integer. Number of initial iterations discarded.
#' @param thin Integer. Thinning interval.
#' @param nChains Integer. Number of independent MCMC chains.
 #' @param var_names Character vector of length \eqn{N} with the names of the endogenous variables.
#' @param lag_names Character vector with labels for lagged regressors (used in summaries/plots).
#' @param country_names Character vector of length \code{length(Y_list)} with group/country names.
#'
#' @return
#' A list with one element per chain, each containing:
#' \itemize{
#'   \item \code{samples} — lists of stored draws for \code{A}, \code{V},
#'   \code{Sigma}, \code{Ag}, \code{Sg}, and scalar hyperparameters;
#'   \item \code{last_state} — the final value of the state matrices;
#'   \item \code{last_hyper} — the final values of scalar hyperparameters.
#' }
#'
#' @details
#' The function internally uses:
#' \itemize{
#'   \item \code{build_starting_values()} for initial states,
#'   \item \code{hyper()} for first hyperparameter draws,
#'   \item \code{step_country_MNIW()} for group-level updates,
#'   \item \code{step_A_MN()}, \code{step_V_IW()}, \code{step_Sigma_W()}
#'         for global parameters,
#'   \item \code{step_hyper_all()} for all scalar hyperparameters.
#' }
#'
#' @examples
#' # Simulate minimal data
#' Y_list <- list(matrix(rnorm(50), 10, 5))
#' Z_list <- list(matrix(rnorm(50), 10, 5))
#'
#' # Variable names (optional)
#' vnames <- paste0("V", 1:5)
#' lnames <- paste0("L", 1:5)
#' cnames <- paste0("C", seq_along(Y_list))
#'
#'
#' out <- run_mcmc(
#'   Y_list, Z_list,
#'   lag = 1,
#'   fixed_list = NULL,
#'   R = 100,
#'   burnin = 50,
#'   thin = 2,
#'   nChains = 1,
#'   var_names = vnames,
#'   lag_names = lnames,
#'   country_names = cnames
#' )
#'
#' # Inspect results
#' str(out)
#'
#' @export
run_mcmc <- function(Y_list, Z_list, lag, fixed_list,
                     R, burnin, thin, nChains,
                     var_names = NULL, lag_names = NULL, country_names = NULL) {

    G <- length(Y_list)
    N <- ncol(Y_list[[1]])
    K <- ncol(Z_list[[1]])

    mcmc <- list(CHAINS = vector("list", nChains))

    STATES  <- vector("list", nChains)
    HYPERS  <- vector("list", nChains)
    SAMPLES <- vector("list", nChains)

    SV    <- build_starting_values(N, K, lag, fixed_list)
    baseSeed <- 123456
    nKeep <- floor(R / thin)
    Rtot  <- R + burnin

    show_pb <- isTRUE(getOption("mcmc.pb", TRUE)) && interactive()
    total_steps <- nChains * Rtot
    if (show_pb) {
        pb <- utils::txtProgressBar(min = 0, max = total_steps, style = 3)
        .step <- 0L
        on.exit(try(close(pb), silent = TRUE), add = TRUE)
    }

    for (c in seq_len(nChains)) {
        set.seed(baseSeed + c)


        HYPERS[[c]] <- hyper(SV)
        STATES[[c]] <- init(SV, Y_list, Z_list, HYPERS[[c]])

        SAMPLES[[c]] <- list(
            A       = vector("list", nKeep),
            V       = vector("list", nKeep),
            Sigma   = vector("list", nKeep),
            Ag      = matrix(vector("list", nKeep * G), nrow = nKeep, ncol = G),
            Sg      = matrix(vector("list", nKeep * G), nrow = nKeep, ncol = G),
            m       = numeric(nKeep),
            s       = numeric(nKeep),
            s_Sigma = numeric(nKeep),
            w       = numeric(nKeep)
        )

        keepIdx <- 0L


        for (i in seq_len(Rtot)) {

            for (g in seq_len(G)) {
                Yg <- Y_list[[g]]
                Zg <- Z_list[[g]]
                STATES[[c]] <- step_country_MNIW(Yg, Zg, SV, STATES[[c]], g)
            }

            STATES[[c]] <- step_A_MN   (SV, STATES[[c]], HYPERS[[c]])
            STATES[[c]] <- step_V_IW   (SV, STATES[[c]], HYPERS[[c]])
            STATES[[c]] <- step_Sigma_W(SV, STATES[[c]])

            HYPERS[[c]] <- step_hyper_all(SV, STATES[[c]], HYPERS[[c]])

            if (i > burnin && ((i - burnin) %% thin == 0L)) {
                keepIdx <- keepIdx + 1L

                for (g in seq_len(G)) {
                    SAMPLES[[c]]$Ag[[keepIdx, g]] <- STATES[[c]]$Ag[[g]]   # KxN
                    SAMPLES[[c]]$Sg[[keepIdx, g]] <- STATES[[c]]$Sg[[g]]   # NxN
                }

                SAMPLES[[c]]$A[[keepIdx]]     <- STATES[[c]]$A           # KxN
                SAMPLES[[c]]$V[[keepIdx]]     <- STATES[[c]]$V           # KxK
                SAMPLES[[c]]$Sigma[[keepIdx]] <- STATES[[c]]$Sigma       # NxN

                SAMPLES[[c]]$m[keepIdx]        <- HYPERS[[c]]$m
                SAMPLES[[c]]$s[keepIdx]        <- HYPERS[[c]]$s
                SAMPLES[[c]]$w[keepIdx]        <- HYPERS[[c]]$w
                SAMPLES[[c]]$s_Sigma[keepIdx]  <- HYPERS[[c]]$s_Sigma
            }

            if (show_pb) {
                .step <- .step + 1L
                if ((.step %% 5L) == 0L || .step == total_steps) {
                    utils::setTxtProgressBar(pb, .step)
                    flush.console()
                }
            }
        }

        chain_out <- list(samples = SAMPLES[[c]],
                          last_state = STATES[[c]],
                          last_hyper = HYPERS[[c]])
        mcmc$CHAINS[[c]] <- chain_out
    }

    if (!is.null(var_names) && !is.null(lag_names) && !is.null(country_names)) {
        mcmc$CHAINS <- rename_output(
            chains        = mcmc$CHAINS,
            var_names     = var_names,
            lag_names     = lag_names,
            country_names = country_names,
            n_groups      = length(country_names),
            R             = R,
            thin          = thin
        )
    }

    return(mcmc)
}
