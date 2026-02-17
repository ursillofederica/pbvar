#' Plot global impulse response functions
#'
#' @description
#' Produces an N x N grid of impulse response plots for the **global**
#' VAR, using posterior draws of IRFs.
#' Each subplot displays the median and equal-tailed 95% credible band.
#'
#' Rows correspond to **responses**,
#' columns correspond to **shocks**.
#'
#' @param Phi_global A 4D array of IRFs:
#'   \code{Phi_global[draw, horizon, i, j]}.
#' @param var_names Character vector of variable names (length N).
#'
#' @export
plot_global_irf <- function(Phi_global, var_names) {

    N <- dim(Phi_global)[3]
    H <- dim(Phi_global)[2]
    horizons <- 0:(H-1)

    par(mfrow=c(N, N), mar=c(3,3,3,1))

    for (i in 1:N) {
        for (j in 1:N) {

            draws_mat <- Phi_global[ , , i, j]
            med  <- apply(draws_mat, 2, median)
            low  <- apply(draws_mat, 2, quantile, probs = 0.025)
            high <- apply(draws_mat, 2, quantile, probs = 0.975)

            plot(
                horizons, med, type="l", lwd=2, col="black",
                ylim=range(c(low, high)),
                xlab="", ylab="",
                main = paste0(var_names[i], " -> shock(", var_names[j], ")")
            )

            polygon(
                c(horizons, rev(horizons)),
                c(low, rev(high)),
                col=adjustcolor("gray80", alpha.f = 0.7),
                border=NA
            )

            lines(horizons, med, lwd=2)
            abline(h=0)
        }
    }

}

