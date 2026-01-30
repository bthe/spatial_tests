# RTMB objective function translating the TMB C++ model logic.

rtmb_Q_spde <- function(spde, kappa) {
  kappa2 <- kappa * kappa
  kappa4 <- kappa2 * kappa2
  kappa4 * spde$M0 + 2 * kappa2 * spde$M1 + spde$M2
}

rtmb_dar1 <- function(x, rho, log = TRUE) {
  if (!log) {
    stop("rtmb_dar1 supports log=TRUE only")
  }
  x <- as.numeric(x)
  rho <- as.numeric(rho)
  if (length(rho) == 0) {
    stop("rtmb_dar1: rho must be numeric")
  }
  n <- length(x)
  if (n == 0) return(0)
  if (n == 1) return(RTMB::dnorm(x, 0, 1, log = TRUE))
  innov_sd <- sqrt(1 - rho * rho)
  RTMB::dnorm(x[1], 0, 1, log = TRUE) +
    sum(RTMB::dnorm(x[-1], rho * x[-n], innov_sd, log = TRUE))
}

rtmb_objective <- function(parms, data) {
  RTMB::getAll(data, parms, warn = FALSE)

  if (!exists("keep", inherits = FALSE)) {
    keep <- rep(1, length(obsVector))
  }

  sigma <- exp(log_sigma)
  kappa <- exp(log_kappa)
  size <- exp(logSize)
  lambda <- exp(log_lambda)

  rho_t <- tan_rho_t
  rho_t[1] <- 2 / (1 + exp(-2 * tan_rho_t[1])) - 1
  rho_l <- tan_rho_l
  rho_l[1] <- 2 / (1 + exp(-2 * tan_rho_l[1])) - 1
  rho_l[2] <- 2 / (1 + exp(-2 * tan_rho_l[2])) - 1
  rho_l[3] <- 2 / (1 + exp(-2 * tan_rho_l[3])) - 1

  if (length(delta_z) > 1) {
    delta_z[-1] <- exp(delta_z[-1])
  }

  Q_s <- rtmb_Q_spde(spdeMatricesS, kappa[1])
  Q_st <- rtmb_Q_spde(spdeMatricesST, kappa[2])

  nll <- 0

  d <- 2
  R <- -log(pcPriorsRange[2]) * pcPriorsRange[1]^(d / 2)
  S <- -log(pcPriorsSD[2]) / pcPriorsSD[1]

  if (spatial == 1) {
    if (usePCpriors == 1) {
      rhoP <- sqrt(8) / kappa[1]
      nll <- nll - log(d / 2 * R * S * rhoP^(-1 - d / 2) *
        exp(-R * rhoP^(-d / 2) - S * sigma[1]))
    }
    f_spatial <- function(x) RTMB::dgmrf(x, 0, Q_s, log = TRUE)
    f_length <- function(x) rtmb_dar1(x, rho_l[1], log = TRUE)
    nll <- nll - RTMB::dseparable(f_spatial, f_length)(xS, log = TRUE)
    scaleS <- 1 / ((4 * pi) * kappa[1] * kappa[1])
    xS <- xS / sqrt(scaleS)
  }

  if (spatioTemporal == 1) {
    if (usePCpriors == 1) {
      rhoP <- sqrt(8) / kappa[2]
      nll <- nll - log(d / 2 * R * S * rhoP^(-1 - d / 2) *
        exp(-R * rhoP^(-d / 2) - S * sigma[2]))
    }
    f_spatial <- function(x) RTMB::dgmrf(x, 0, Q_st, log = TRUE)
    f_time <- function(x) rtmb_dar1(x, rho_t[1], log = TRUE)
    f_length <- function(x) rtmb_dar1(x, rho_l[2], log = TRUE)
    nll <- nll - RTMB::dseparable(f_spatial, f_time, f_length)(xST, log = TRUE)
    scaleST <- 1 / ((4 * pi) * kappa[2] * kappa[2])
    xST <- xST / sqrt(scaleST)
  }

  if (useNugget == 1) {
    nHaul <- as.integer(sum(nStationsEachYear))
    Q_nuggetIID <- Matrix::Diagonal(nHaul)
    f_haul <- function(x) RTMB::dgmrf(x, 0, Q_nuggetIID, log = TRUE)
    f_length <- function(x) rtmb_dar1(x, rho_l[3], log = TRUE)
    nll <- nll - RTMB::dseparable(f_length, f_haul)(nugget, log = TRUE)

    Q_nuggetIIDI <- Matrix::Diagonal(length(xInt))
    f_int <- function(x) RTMB::dgmrf(x, 0, Q_nuggetIIDI, log = TRUE)
    nll <- nll - RTMB::dseparable(f_length, f_int)(nuggetIndex, log = TRUE)
  }

  nSplineDepth <- length(betaDepth) / 2
  parDepth1 <- betaDepth[seq_len(nSplineDepth)]
  parDepth2 <- betaDepth[seq_len(nSplineDepth) + nSplineDepth]
  if (splineDepth == 1) {
    nll <- nll - 0.5 * Sdim * log_lambda[1] - 0.5 * lambda[1] * as.numeric(t(parDepth1) %*% (S_depth %*% parDepth1))
  } else if (splineDepth == 2) {
    nll <- nll - 0.5 * Sdim * log_lambda[1] - 0.5 * lambda[1] * as.numeric(t(parDepth1) %*% (S_depth %*% parDepth1))
    nll <- nll - 0.5 * Sdim * log_lambda[2] - 0.5 * lambda[2] * as.numeric(t(parDepth2) %*% (S_depth %*% parDepth2))
  }

  betaSunLow <- betaSun[seq_len(nBasisSunAlt * 2)]
  betaSunHigh <- betaSun[seq_len(nBasisSunAlt * 2) + nBasisSunAlt * 2]

  maxStations <- max(nStationsEachYear)
  mu <- matrix(0, nrow(fishObsMatrix), ncol(fishObsMatrix))
  deltaMatrixS <- matrix(0, nrow = numberOfLengthGroups, ncol = maxStations)
  deltaMatrixST <- matrix(0, nrow = numberOfLengthGroups, ncol = maxStations)

  timeInDayLow <- X_sunAlt %*% betaSunLow
  timeInDayHigh <- X_sunAlt %*% betaSunHigh
  depthEffect1 <- X_depth %*% parDepth1
  depthEffect2 <- X_depth %*% parDepth2

  validation <- rep(0, numberOfLengthGroups)
  counter <- 1
  for (y in seq_len(length(nStationsEachYear))) {
    As <- A_ListS[[y]]
    Ast <- A_ListST[[y]]

    if (as.integer(lengthGroupsReduced[1]) == as.integer(lengthGroupsReduced[2])) {
      for (l in seq_len(numberOfLengthGroups)) {
        deltaS <- As %*% xS[, as.integer(lengthGroupsReduced[l])]
        deltaS2 <- As %*% xS[, as.integer(lengthGroupsReduced[l]) + 1]
        deltaST <- Ast %*% xST[, y, as.integer(lengthGroupsReduced[l])]
        deltaST2 <- Ast %*% xST[, y, as.integer(lengthGroupsReduced[l]) + 1]

        for (s in seq_len(nStationsEachYear[y])) {
          deltaMatrixS[l, s] <- weigthLength[l] * deltaS[s] + (1 - weigthLength[l]) * deltaS2[s]
          deltaMatrixST[l, s] <- weigthLength[l] * deltaST[s] + (1 - weigthLength[l]) * deltaST2[s]
        }
      }
    } else {
      for (l in seq_len(numberOfLengthGroups)) {
        deltaS <- As %*% xS[, as.integer(lengthGroupsReduced[l])]
        deltaST <- Ast %*% xST[, y, as.integer(lengthGroupsReduced[l])]
        for (s in seq_len(nStationsEachYear[y])) {
          deltaMatrixS[l, s] <- deltaS[s]
          deltaMatrixST[l, s] <- deltaST[s]
        }
      }
    }

    for (s in seq_len(nStationsEachYear[y])) {
      for (l in seq_len(numberOfLengthGroups)) {
        covariatesConvexW <- (numberOfLengthGroups - l) / (numberOfLengthGroups - 1)
        mu[counter, l] <- exp(beta0[y, l] +
          covariatesConvexW * timeInDayLow[counter] + (1 - covariatesConvexW) * timeInDayHigh[counter] +
          deltaMatrixS[l, s] * sigma[1] +
          deltaMatrixST[l, s] * sigma[2] +
          betahaul[l] * log(dist[counter]) +
          covariatesConvexW * depthEffect1[counter] + (1 - covariatesConvexW) * depthEffect2[counter] +
          nugget[l, counter] * sigma[3])

        log_var_minus_mu <- log(mu[counter, l] * mu[counter, l] * size)
        if (predMatrix[counter, l] == 0) {
          if (zeroInflated == 1) {
            muZero <- exp(delta_z[1] +
              delta_z[2] * beta0[y, l] +
              delta_z[3] * (covariatesConvexW * timeInDayLow[counter] + (1 - covariatesConvexW) * timeInDayHigh[counter]) +
              delta_z[4] * deltaMatrixS[l, s] * sigma[1] +
              delta_z[5] * deltaMatrixST[l, s] * sigma[2] +
              delta_z[6] * (covariatesConvexW * depthEffect1[counter] + (1 - covariatesConvexW) * depthEffect2[counter]) +
              delta_z[8] * nugget[l, counter] * sigma[3] +
              delta_z[9] * betahaul[l] * log(dist[counter]))
            pZero <- dpois(0, muZero, log = TRUE)
            if (obsModel == 1) {
              pPos <- RTMB::dnbinom_robust(obsVector[idxStart[counter] + l + 1], log(mu[counter, l]), log_var_minus_mu, log = TRUE) +
                RTMB::logspace_sub(0, pZero)
            } else if (obsModel == 2) {
              pPos <- dpois(obsVector[idxStart[counter] + l + 1], mu[counter, l], log = TRUE) +
                RTMB::logspace_sub(0, pZero)
            } else {
              stop("obsModel not implemented")
            }

            if (fishObsMatrix[counter, l] == 0) {
              nll <- nll - keep[idxStart[counter] + l + 1] * RTMB::logspace_add(pZero, pPos)
            } else {
              nll <- nll - keep[idxStart[counter] + l + 1] * pPos
            }
          } else {
            if (obsModel == 1) {
              nll <- nll - keep[idxStart[counter] + l + 1] *
                RTMB::dnbinom_robust(obsVector[idxStart[counter] + l + 1], log(mu[counter, l]), log_var_minus_mu, log = TRUE)
            } else if (obsModel == 2) {
              nll <- nll - keep[idxStart[counter] + l + 1] *
                dpois(obsVector[idxStart[counter] + l + 1], mu[counter, l], log = TRUE)
            }
          }
        } else {
          validation[l] <- validation[l] + mu[counter, l]
        }
      }
      counter <- counter + 1
    }
  }

  if (doValidation == 1) {
    RTMB::ADREPORT(validation)
  }

  nYears <- length(nStationsEachYear)
  muReport <- array(0, dim = c(nYears, nStrata, ncol(fishObsMatrix)))
  muReportFull <- matrix(0, nYears, ncol(fishObsMatrix))
  muReportSelected <- matrix(0, nYears, ncol(fishObsMatrix))
  COGx <- matrix(0, nYears, ncol(fishObsMatrix))
  COGy <- matrix(0, nYears, ncol(fishObsMatrix))
  muReportSelectedExp <- matrix(0, nYears, ncol(fishObsMatrix))
  muReportFullExp <- matrix(0, nYears, ncol(fishObsMatrix))

  depthEffectInt1 <- X_depth_int %*% parDepth1
  depthEffectInt2 <- X_depth_int %*% parDepth2
  timeInDayEffectIntLow <- X_sunAltIntegrate %*% betaSunLow
  timeInDayEffectIntHigh <- X_sunAltIntegrate %*% betaSunHigh

  for (y in seq_len(nYears)) {
    for (l in seq_len(numberOfLengthGroups)) {
      covariatesConvexW <- (numberOfLengthGroups - l) / (numberOfLengthGroups - 1)
      if (as.integer(lengthGroupsReduced[1]) == as.integer(lengthGroupsReduced[2])) {
        deltaPredS <- weigthLength[l] * (ApredS %*% xS[, as.integer(lengthGroupsReduced[l])]) +
          (1 - weigthLength[l]) * (ApredS %*% xS[, as.integer(lengthGroupsReduced[l]) + 1])
        deltaPredST <- weigthLength[l] * (ApredST %*% xST[, y, as.integer(lengthGroupsReduced[l])]) +
          (1 - weigthLength[l]) * (ApredST %*% xST[, y, as.integer(lengthGroupsReduced[l]) + 1])
      } else {
        deltaPredS <- ApredS %*% xS[, as.integer(lengthGroupsReduced[l])]
        deltaPredST <- ApredST %*% xST[, y, as.integer(lengthGroupsReduced[l])]
      }

      COGtotal <- 0
      for (strata in seq_len(nStrata)) {
        predObsInStrata <- predLoc[strata, ]
        idx <- predObsInStrata[predObsInStrata >= 0] + 1
        muReport[y, strata, l] <- 0
        if (length(idx) > 0) {
          for (i in idx) {
            muThis <- exp(beta0[y, l] +
              covariatesConvexW * timeInDayEffectIntLow[1] + (1 - covariatesConvexW) * timeInDayEffectIntHigh[1] +
              covariatesConvexW * depthEffectInt1[i] + (1 - covariatesConvexW) * depthEffectInt2[i] +
              deltaPredS[i] * sigma[1] +
              deltaPredST[i] * sigma[2] +
              nuggetIndex[l, i] * sigma[3])

            muReport[y, strata, l] <- muReport[y, strata, l] + muThis
            COGx[y, l] <- COGx[y, l] + muThis * xInt[i]
            COGy[y, l] <- COGy[y, l] + muThis * yInt[i]
            COGtotal <- COGtotal + muThis
          }
          nIntegrate <- length(idx)
          muReportFull[y, l] <- muReportFull[y, l] + muReport[y, strata, l] * areas[strata] / nIntegrate
          if (any(as.integer(selectedStratas) == strata)) {
            muReportSelected[y, l] <- muReportSelected[y, l] + muReport[y, strata, l] * areas[strata] / nIntegrate
          }
          muReport[y, strata, l] <- log(muReport[y, strata, l] * areas[strata] / nIntegrate)
        }
      }

      muReportFullExp[y, l] <- muReportFull[y, l]
      muReportSelectedExp[y, l] <- muReportSelected[y, l]
      muReportFull[y, l] <- log(muReportFull[y, l])
      muReportSelected[y, l] <- log(muReportSelected[y, l])

      COGx[y, l] <- COGx[y, l] / COGtotal
      COGy[y, l] <- COGy[y, l] / COGtotal
    }
  }

  if (doDetailedRep == 1) {
    RTMB::ADREPORT(muReportFullExp)
    RTMB::ADREPORT(muReportSelectedExp)
    RTMB::ADREPORT(muReportFull)
    RTMB::ADREPORT(muReportSelected)
    RTMB::ADREPORT(COGx)
    RTMB::ADREPORT(COGy)
  } else if (doDetailedRep == 2) {
    RTMB::ADREPORT(muReportFullExp)
    RTMB::ADREPORT(muReportSelectedExp)
    RTMB::ADREPORT(muReportFull)
    RTMB::ADREPORT(muReportSelected)
    RTMB::ADREPORT(COGx)
    RTMB::ADREPORT(COGy)

    fourierReportLow <- X_sunAltReport %*% betaSunLow
    fourierReportHigh <- X_sunAltReport %*% betaSunHigh
    RTMB::ADREPORT(fourierReportLow)
    RTMB::ADREPORT(fourierReportHigh)

    depthReport1 <- X_depthReport %*% parDepth1
    depthReport2 <- X_depthReport %*% parDepth2
    RTMB::ADREPORT(depthReport1)
    RTMB::ADREPORT(depthReport2)
  }

  nll
}
