#' Summarize an \code{mcmc.list} object (internal)
#'
#' @description
#' Computes basic posterior summary statistics
#' stored in an \code{mcmc.list} object.
#' The function extracts:
#' \itemize{
#'   \item posterior mean and standard deviation;
#'   \item central credible intervals (2.5%, 50%, 97.5%);
#'   \item 95% HPD interval;
#'   \item Gelman–Rubin R-hat diagnostic;
#'   \item effective sample size (ESS).
#' }
#'
#' This function is internal and used by \code{summary_par()}.
#'
#' @param mcmc_list An object of class \code{mcmc.list} (from the \pkg{coda} package)
#'   containing samples of a single scalar parameter across chains.
#'
#' @return A one-row \code{data.frame} with posterior summary statistics.
#'
#' @keywords internal
#' @noRd
summarize_mcmc <- function(mcmc_list) {

    sm <- summary(mcmc_list)

    mean  <- sm$statistics["Mean"]
    sd    <- sm$statistics["SD"]
    q025  <- sm$quantiles["2.5%"]
    q50   <- sm$quantiles["50%"]
    q975  <- sm$quantiles["97.5%"]

    hpd      <- coda::HPDinterval(mcmc_list, prob = 0.95)
    HPD_low  <- hpd[[1]][1]
    HPD_high <- hpd[[1]][2]

    Rhat <- coda::gelman.diag(mcmc_list, autoburnin=FALSE)$psrf[1]
    ESS  <- coda::effectiveSize(mcmc_list)

    data.frame(
        mean = mean,
        sd = sd,
        q025 = q025,
        q50 = q50,
        q975 = q975,
        HPD_low = HPD_low,
        HPD_high = HPD_high,
        Rhat = Rhat,
        ESS = ESS
    )
}

#' Convert a list of summary data frames into a tidy table (internal)
#'
#' @description
#' Takes a named list of posterior summary data frames and combines them
#' into a single tidy data frame.
#'
#' This helper is used by \code{summary_par()}.
#'
#' @param summary_list A named list where each element is a one-row
#'   \code{data.frame} produced by \code{summarize_mcmc()}.
#'
#' @return A tidy \code{data.frame} with one row per parameter and columns
#'   containing posterior summary statistics.
#'
#' @keywords internal
#' @noRd
tidy_summary <- function(summary_list) {

    out <- lapply(names(summary_list), function(name) {
        df <- summary_list[[name]]
        df$param <- name
        df
    })

    out <- do.call(rbind, out)
    out <- out[, c("param", setdiff(names(out), "param"))]

    rownames(out) <- NULL
    out
}
