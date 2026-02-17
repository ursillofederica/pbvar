# pbvar

Small R package developed for a university project.
It implements a Bayesian hierarchical Panel VAR (PBVAR) estimated via MCMC, with utilities for posterior summaries and impulse response functions.

## What’s inside

- Hierarchical PBVAR sampler
- Posterior summaries and basic diagnostics
- Global impulse response functions
- Example dataset (G7 panel)

## Install

```r
# install.packages("remotes")
remotes::install_github("ursillofederica/pbvar")
```

## Quick example

```r
library(pbvar)

data(panel_g7_bvar)

# Run a small MCMC (for real results increase iterations and nChains)
out <- run_mcmc(
  Y_list = panel_g7_bvar$Y_list,
  Z_list = panel_g7_bvar$Z_list,
  lag = panel_g7_bvar$lag,
  fixed_list = panel_g7_bvar$fixed_list,
  R = 500,
  burnin = 200,
  thin = 2,
  nChains = 2,
  var_names = panel_g7_bvar$var_names,
  lag_names = panel_g7_bvar$lag_names,
  country_names = panel_g7_bvar$country_names
)

# Posterior summaries
par <- summary_par(
  out$CHAINS,
  var_names = panel_g7_bvar$var_names,
  lag_names = panel_g7_bvar$lag_names,
  country_names = panel_g7_bvar$country_names,
  par = "all"
)

# Global IRF
irf <- compute_global_irf(out)
plot_global_irf(irf)
```

## Vignette

See `vignettes/example_G7.Rmd` for a full example.

## Notes

This code was written for learning purposes and is not meant to be a polished CRAN package.
