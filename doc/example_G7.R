## -----------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  out.width = "90%"
)

## ----setup--------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(pbvar)
if (!requireNamespace("coda", quietly = TRUE)) {
  stop("Package 'coda' is required to build this vignette. Install it with install.packages('coda').")
}
if (!requireNamespace("lattice", quietly = TRUE)) {
  stop("Package 'lattice' is required to build this vignette. Install it with install.packages('lattice').")
}

## -----------------------------------------------------------------------------
data(panel_g7_bvar)

# USA
# Select the raw data for USA
df_USA <- panel_g7_bvar$USA

# Transformations:
# GDP_USA:   log-difference × 100   (quarterly growth rate)
# CPI_USA:   first difference       (already expressed in rates)
# U_USA:     first difference       (changes in unemployment)
df_USA$GDPd <- c(NA, diff(log(df_USA$GDP_USA)) * 100)
df_USA$CPId <- c(NA, diff(df_USA$CPI_USA))
df_USA$Ud   <- c(NA, diff(df_USA$U_USA))

# Remove first two rows (contain NA created by differencing/log-diff)
df_USA <- df_USA[-c(1, 2), ]

# Construct matrix of transformed variables in VAR order
# Here: GDPd (growth), CPId (inflation), Ud (unemployment changes)
M_USA <- as.matrix(df_USA[, c("GDPd", "CPId", "Ud")])
T_USA <- nrow(M_USA)

# Build Y_g and Z_g for a VAR(1):
# Y_g is X_{t}, Z_g is X_{t-1}
Y_USA <- M_USA[2:T_USA, ]
X_USA <- M_USA[1:(T_USA - 1), ]

# CAN
df_CAN <- panel_g7_bvar$CAN

df_CAN$GDPd <- c(NA, diff(log(df_CAN$GDP_CAN)) * 100)
df_CAN$CPId <- c(NA, diff(df_CAN$CPI_CAN))
df_CAN$Ud   <- c(NA, diff(df_CAN$U_CAN))

df_CAN <- df_CAN[-c(1, 2), ]

M_CAN <- as.matrix(df_CAN[, c("GDPd", "CPId", "Ud")])
T_CAN <- nrow(M_CAN)

Y_CAN <- M_CAN[2:T_CAN, ]
X_CAN <- M_CAN[1:(T_CAN - 1), ]

# GBR
df_GBR <- panel_g7_bvar$GBR

df_GBR$GDPd <- c(NA, diff(log(df_GBR$GDP_GBR)) * 100)
df_GBR$CPId <- c(NA, diff(df_GBR$CPI_GBR))
df_GBR$Ud   <- c(NA, diff(df_GBR$U_GBR))

df_GBR <- df_GBR[-c(1, 2), ]

M_GBR <- as.matrix(df_GBR[, c("GDPd", "CPId", "Ud")])
T_GBR <- nrow(M_GBR)

Y_GBR <- M_GBR[2:T_GBR, ]
X_GBR <- M_GBR[1:(T_GBR - 1), ]

# DEU
df_DEU <- panel_g7_bvar$DEU

df_DEU$GDPd <- c(NA, diff(log(df_DEU$GDP_DEU)) * 100)
df_DEU$CPId <- c(NA, diff(df_DEU$CPI_DEU))
df_DEU$Ud   <- c(NA, diff(df_DEU$U_DEU))

df_DEU <- df_DEU[-c(1, 2), ]

M_DEU <- as.matrix(df_DEU[, c("GDPd", "CPId", "Ud")])
T_DEU <- nrow(M_DEU)

Y_DEU <- M_DEU[2:T_DEU, ]
X_DEU <- M_DEU[1:(T_DEU - 1), ]

# FRA
df_FRA <- panel_g7_bvar$FRA

df_FRA$GDPd <- c(NA, diff(log(df_FRA$GDP_FRA)) * 100)
df_FRA$CPId <- c(NA, diff(df_FRA$CPI_FRA))
df_FRA$Ud   <- c(NA, diff(df_FRA$U_FRA))

df_FRA <- df_FRA[-c(1, 2), ]

M_FRA <- as.matrix(df_FRA[, c("GDPd", "CPId", "Ud")])
T_FRA <- nrow(M_FRA)

Y_FRA <- M_FRA[2:T_FRA, ]
X_FRA <- M_FRA[1:(T_FRA - 1), ]

# ITA
df_ITA <- panel_g7_bvar$ITA

df_ITA$GDPd <- c(NA, diff(log(df_ITA$GDP_ITA)) * 100)
df_ITA$CPId <- c(NA, diff(df_ITA$CPI_ITA))
df_ITA$Ud   <- c(NA, diff(df_ITA$U_ITA))

df_ITA <- df_ITA[-c(1, 2), ]

M_ITA <- as.matrix(df_ITA[, c("GDPd", "CPId", "Ud")])
T_ITA <- nrow(M_ITA)

Y_ITA <- M_ITA[2:T_ITA, ]
X_ITA <- M_ITA[1:(T_ITA - 1), ]

# JPN
df_JPN <- panel_g7_bvar$JPN

df_JPN$GDPd <- c(NA, diff(log(df_JPN$GDP_JPN)) * 100)
df_JPN$CPId <- c(NA, diff(df_JPN$CPI_JPN))
df_JPN$Ud   <- c(NA, diff(df_JPN$U_JPN))

df_JPN <- df_JPN[-c(1, 2), ]

M_JPN <- as.matrix(df_JPN[, c("GDPd", "CPId", "Ud")])
T_JPN <- nrow(M_JPN)

Y_JPN <- M_JPN[2:T_JPN, ]
X_JPN <- M_JPN[1:(T_JPN - 1), ]

Yg_list <- list(
  Y_USA,Y_CAN,Y_GBR,Y_DEU,
  Y_FRA,Y_ITA,Y_JPN
)

Xg_list <- list(
  X_USA,X_CAN,X_GBR,X_DEU,
  X_FRA,X_ITA,X_JPN
)


## -----------------------------------------------------------------------------
Yg_list <- unname(Yg_list)
Xg_list <- unname(Xg_list)

# Demean each country-specific block
Yg_list_demean <- lapply(Yg_list, function(Y) scale(Y, center = TRUE, scale = FALSE))
Xg_list_demean <- lapply(Xg_list, function(X) scale(X, center = TRUE, scale = FALSE))

## -----------------------------------------------------------------------------
lag        <- 1
fixed_list <- NULL
R          <- 2000
burnin     <- 500
thin       <- 5
nChains    <- 2

var_names     <- c("GDP", "CPI", "U")
lag_names     <- paste0(var_names, "_lag1")
country_names <- c("USA","CAN","GBR","DEU","FRA","ITA","JPN")

mcmc_out <- run_mcmc(Yg_list_demean, Xg_list_demean, lag, fixed_list,
                     R, burnin, thin, nChains,
                     var_names, lag_names, country_names)

## -----------------------------------------------------------------------------
chains  <- mcmc_out$CHAINS
nChains <- length(chains)

var_names     <- c("GDP", "CPI", "U")
lag_names     <- paste0(var_names, "_lag1")
country_names <- c("USA","CAN","GBR","DEU","FRA","ITA","JPN")

par <- summary_par(chains,var_names,lag_names,country_names, par = "all")

tabA    <-par$A
tabAg   <-par$Ag
tabSg   <-par$Sg
tabSigma<-par$Sigma
tabV    <-par$V


## -----------------------------------------------------------------------------
# A
mean(abs(tabA$Rhat-1))
sum(abs(tabA$Rhat-1)>0.1)

#Ag
mean(abs(tabAg$Rhat-1))
sum(abs(tabAg$Rhat-1)>0.1)

#Sg
mean(abs(tabSg$Rhat-1))
sum(abs(tabSg$Rhat-1)>0.1)

#Sigma
mean(abs(tabSigma$Rhat-1))
sum(abs(tabSigma$Rhat-1)>0.1)

#V
mean(abs(tabV$Rhat-1))
sum(abs(tabV$Rhat-1)>0.1)

## -----------------------------------------------------------------------------
tot_ESS<-c(tabA$ESS,tabAg$ESS,tabSg$ESS,
           tabSigma$ESS,tabV$ESS)
plot(tot_ESS,type='l')
abline(v=9,col='red')
abline(v=(9+63),col='purple')
abline(v=(9+63+63),col='dark blue')
abline(v=(9+9+63+63),col='dark orange')
text(3,1200,'A')
text(42,1200,'Ag')
text(100,1200,'Sg')
text(139.5,1200,'Sigma')
text(150,1200,'V')
abline(h=min(tot_ESS))
text(which.min(tot_ESS)+1, min(tot_ESS)-.5,'min: 292.6293',cex=.7)

## -----------------------------------------------------------------------------
# Positive definite check 
mat<-list(
  meanV<-matrix(tabV[,2],ncol=3),
  meanSig<-matrix(tabSigma[,2],ncol=3),
  meanSg<-array(tabSg[,2],c(3,3,7)))

defpos<-function(mat,indx){
  if(is.na(dim(mat[[indx]])[3])){
    if(sum(eigen(mat[[indx]])$values>0)!=3){
      cat("There is a problem in the matrix", indx,'\n',
         'The eigenvalues of the matrix are', eigen(mat[[indx]])$values,'\n')}
    else cat('The matrix is positive definite')
  }
  else{
    for(j in 1:7){
      if(sum(eigen(mat[[indx]][,,j])$values>0)!=3){
        cat("There is a problem in the matrix", i,'\n',
            'The eigenvalues of the group matrix',j,' are', eigen(mat[[indx]][,,j])$values,'\n')
      }
      else cat('The matrix is positive definite \n')
    }
  }
  
}

defpos(mat,1)
defpos(mat,2)
defpos(mat,3)

# unit roots check
K <- 3  
G <- length(country_names)

A_mean <- matrix(tabA$mean, nrow = K, byrow = TRUE)
A_mean

ev_A   <- eigen(A_mean)$values
mod_A  <- Mod(ev_A)
ev_A
mod_A
max(mod_A)  

max_mod_g <- c(G)
for (g in seq_along(country_names)) {

  idx <- ((g - 1) * K * K + 1):(g * K * K) 
  Ag_mean_g <- matrix(tabAg$mean[idx], nrow = K, byrow = TRUE)

  ev  <- eigen(Ag_mean_g)$values
  mod <- Mod(ev)
  
  max_mod_g[g] <- max(mod)
}

max_mod_g

## -----------------------------------------------------------------------------

tot_Ag <- array(tabAg$mean,dim=c(3,3,7))
matA   <-matrix(tabA$mean,ncol=3)

colnames(matA) <- c('GDP_lag1','CPI_lag1','U_lag1')
rownames(matA) <- c('GDP','CPI','U')

dimnames(tot_Ag) <- list(
  rownames = c('GDP','CPI','U'),
  colnames = c('GDP_lag1','CPI_lag1','U_lag1'),
  znames   = c("USA", "CAN","GBR","DEU","FRA","ITA","JPN")
)

apply(tot_Ag,3,function(x) sqrt(sum((x-matA)^2)))

## ----fig.width=6, fig.height=6------------------------------------------------
H <- 10

global_irfs <- compute_global_irf(chains,H)

plot_global_irf(global_irfs, var_names)

