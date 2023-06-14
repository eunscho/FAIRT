study6 <- function(dich = T, start_rep = 1, end_rep = 20) {
  library(tidyverse)
  library(mirt)
  numcat <- ifelse(dich, 2, 5)
  models <- 1:10
  stats <- c("MAE", "RMSE")
  out_ncol <- length(stats) * length(models) + 3 # 3: numcat, rep_set, rep
  optrmse <- ifelse(dich, 1.739, 1.753)
  optmae <- ifelse(dich, 1.734, 1.746)
  rep_per_sets <- 1:1000
  name <- "study6_"
  rep_sets <- start_rep:end_rep
  #========================================================================
  # Generate
  #========================================================================
  generate <- function(numcat, model, rep_set, rep) {
    # seed and variable names
    set.seed(10000 * numcat + 1000 * rep_set + rep)
    NITEM <- 10
    n <- 10000
    skew <- 0
    kurt <- 0
    theta <- semTools::mvrnonnorm(n, c(0,0), matrix(c(1, 0, 0, 1), 2, 2),
                                  skewness = c(skew, 0),
                                  kurtosis = c(kurt, 0))[, 1]
    alpha <- runif(NITEM, .6, 1.4)
    lambda <- alpha / sqrt(1 + alpha ^ 2)
    psi <- sqrt(1 - lambda^2)
    score <-  matrix(vector("double", n * NITEM), n, NITEM)
    if (model == 1) { # FA model
      load <- alpha / sqrt(1 + alpha ^ 2)
      e <- matrix(rnorm(n * NITEM), nrow = n, ncol = NITEM)
      score_con <- theta %*% t(load) + e %*% sqrt(diag(NITEM) - diag(diag(load %*% t(load))))
    } else if (model == 3) {
      D <- 1.702
    } else if (model == 7) {
      D <- ifelse(numcat == 2, 1.734, 1.74)
    } else if (model == 8) {
      D <- ifelse (numcat == 2, 1.739, 1.746)
    } else if (model == 9) {
      D <- ifelse (numcat == 2, 1.75, 1.753)
    } else if (model == 10) {
      D <- 1.76
    } else if (model != 2) { # 4, 5, 6
      D <- 1.67 + .01 * model
      # D <- 1.66 + .01 * model
    }
    if (numcat == 2) {
      beta <- runif(NITEM, -1.3, 1.3)
      if (model == 1) { # FA model
        thresh <- beta / sqrt(1 + alpha ^ 2)
        for (i in 1:n) {
          for (j in 1:NITEM) {
            score[i, j] <- ifelse(score_con[i, j] < thresh[j], 0, 1)
          }
        }
      } else { # IRT model
        probs <- matrix(vector("double", n * NITEM), n, NITEM)
        among <- c(0, 1)
        for (i in 1:n) {
          for (j in 1:NITEM) {
            if (model == 2) { # normal ogive model
              probs[i, j] <- pnorm(alpha[j] * theta[i] - beta[j])
            } else { # logistic model
              probs[i, j] <- LaplacesDemon::invlogit(D * alpha[j] * theta[i] - D * beta[j])
            }
            score[i, j] <- sample(among, size = 1, replace = TRUE,
                                  prob = c(1 - probs[i, j], probs[i, j]))
          }
        }
      }
    } else {
      beta <- matrix(vector("double", 4 * NITEM), nrow = NITEM, ncol = 4)
      score <-  matrix(vector("double", n * NITEM), n, NITEM)
      for (i in 1:NITEM) {
        beta[i, ] <- sort(runif(4, -2, 2))
        beta[i, 1] <- beta[i, 1] - .2
        beta[i, 2] <- beta[i, 2] - .1
        beta[i, 3] <- beta[i, 3] + .1
        beta[i, 4] <- beta[i, 4] + .2
      }
      if (model == 1) { # FA model
        thresh <- matrix(vector("double", 4 * NITEM), nrow = NITEM, ncol = 4)
        for (i in 1:NITEM) {
          for (j in 1:4) {
            thresh[i, j] <- beta[i, j] / sqrt(1 + alpha[i] ^ 2)
          }
        }
        for (i in 1:n) {
          for (j in 1:NITEM) {
            score[i, j] <- ifelse(score_con[i, j] < thresh[j, 1], 0,
                                  ifelse(score_con[i, j] < thresh[j, 2], 1,
                                         ifelse(score_con[i, j] < thresh[j, 3], 2,
                                                ifelse(score_con[i, j] < thresh[j, 4], 3, 4))))        }
        }
      } else { # IRT model
        probs <- array(vector("double", n * NITEM * 4), dim = c(n, NITEM, 4))
        among <- c(0, 1, 2, 3, 4)
        for (i in 1:n) {
          for (j in 1:NITEM) {
            for (k in 1:4) {
              if (model == 2) { # FA model
                probs[i, j, k] <- pnorm(alpha[j] * theta[i] - beta[j, k])
              } else { # logistic model
                probs[i, j, k] <- LaplacesDemon::invlogit(D * alpha[j] * theta[i] - D * beta[j, k])
              }
            }
            score[i, j] <- sample(among, size = 1, replace = TRUE,
                                  prob = c(1 - probs[i, j, 1],
                                           probs[i, j, 1] - probs[i, j, 2],
                                           probs[i, j, 2] - probs[i, j, 3],
                                           probs[i, j, 3] - probs[i, j, 4],
                                           probs[i, j, 4]))
          }
        }
      }
    }
    colnames(score) <- paste0("X", 1:10)
    out <- list(dat = score, alpha = alpha, beta = beta)
    return(out)
  }
  #========================================================================
  # Analyze
  #========================================================================
  analyze <- function(numcat, genout) {
    library(mirt)
    dat <- genout$dat
    alpha <- genout$alpha
    beta <- genout$beta
    res <- mirt(dat, model = 1, method = "EM")
    alpha1 <- alpha[1]
    alpha2 <- alpha[2]
    alpha3 <- alpha[3]
    alpha4 <- alpha[4]
    alpha5 <- alpha[5]
    alpha6 <- alpha[6]
    alpha7 <- alpha[7]
    alpha8 <- alpha[8]
    alpha9 <- alpha[9]
    alpha10 <- alpha[10]
    a1 <- coef(res)$X1[1]
    a2 <- coef(res)$X2[1]
    a3 <- coef(res)$X3[1]
    a4 <- coef(res)$X4[1]
    a5 <- coef(res)$X5[1]
    a6 <- coef(res)$X6[1]
    a7 <- coef(res)$X7[1]
    a8 <- coef(res)$X8[1]
    a9 <- coef(res)$X9[1]
    a10 <- coef(res)$X10[1]
    if (numcat == 2) {
      betai1c1 <- beta[1]
      betai2c1 <- beta[2]
      betai3c1 <- beta[3]
      betai4c1 <- beta[4]
      betai5c1 <- beta[5]
      betai6c1 <- beta[6]
      betai7c1 <- beta[7]
      betai8c1 <- beta[8]
      betai9c1 <- beta[9]
      betai10c1 <- beta[10]
      bi1c1 <- coef(res)$X1[2]
      bi2c1 <- coef(res)$X2[2]
      bi3c1 <- coef(res)$X3[2]
      bi4c1 <- coef(res)$X4[2]
      bi5c1 <- coef(res)$X5[2]
      bi6c1 <- coef(res)$X6[2]
      bi7c1 <- coef(res)$X7[2]
      bi8c1 <- coef(res)$X8[2]
      bi9c1 <- coef(res)$X9[2]
      bi10c1 <- coef(res)$X10[2]
      betai1c2 <- betai2c2 <- betai3c2 <- betai4c2 <- betai5c2 <- NA
      betai6c2 <- betai7c2 <- betai8c2 <- betai9c2 <- betai10c2 <- NA
      betai1c3 <- betai2c3 <- betai3c3 <- betai4c3 <- betai5c3 <- NA
      betai6c3 <- betai7c3 <- betai8c3 <- betai9c3 <- betai10c3 <- NA
      betai1c4 <- betai2c4 <- betai3c4 <- betai4c4 <- betai5c4 <- NA
      betai6c4 <- betai7c4 <- betai8c4 <- betai9c4 <- betai10c4 <- NA
      bi1c2 <- bi2c2 <- bi3c2 <- bi4c2 <- bi5c2 <- NA
      bi6c2 <- bi7c2 <- bi8c2 <- bi9c2 <- bi10c2 <- NA
      bi1c3 <- bi2c3 <- bi3c3 <- bi4c3 <- bi5c3 <- NA
      bi6c3 <- bi7c3 <- bi8c3 <- bi9c3 <- bi10c3 <- NA
      bi1c4 <- bi2c4 <- bi3c4 <- bi4c4 <- bi5c4 <- NA
      bi6c4 <- bi7c4 <- bi8c4 <- bi9c4 <- bi10c4 <- NA
    } else {
      betai1c1 <- beta[1, 1]
      betai2c1 <- beta[2, 1]
      betai3c1 <- beta[3, 1]
      betai4c1 <- beta[4, 1]
      betai5c1 <- beta[5, 1]
      betai6c1 <- beta[6, 1]
      betai7c1 <- beta[7, 1]
      betai8c1 <- beta[8, 1]
      betai9c1 <- beta[9, 1]
      betai10c1 <- beta[10, 1]
      betai1c2 <- beta[1, 2]
      betai2c2 <- beta[2, 2]
      betai3c2 <- beta[3, 2]
      betai4c2 <- beta[4, 2]
      betai5c2 <- beta[5, 2]
      betai6c2 <- beta[6, 2]
      betai7c2 <- beta[7, 2]
      betai8c2 <- beta[8, 2]
      betai9c2 <- beta[9, 2]
      betai10c2 <- beta[10, 3]
      betai1c3 <- beta[1, 3]
      betai2c3 <- beta[2, 3]
      betai3c3 <- beta[3, 3]
      betai4c3 <- beta[4, 3]
      betai5c3 <- beta[5, 3]
      betai6c3 <- beta[6, 3]
      betai7c3 <- beta[7, 3]
      betai8c3 <- beta[8, 3]
      betai9c3 <- beta[9, 3]
      betai10c3 <- beta[10, 3]
      betai1c4 <- beta[1, 4]
      betai2c4 <- beta[2, 4]
      betai3c4 <- beta[3, 4]
      betai4c4 <- beta[4, 4]
      betai5c4 <- beta[5, 4]
      betai6c4 <- beta[6, 4]
      betai7c4 <- beta[7, 4]
      betai8c4 <- beta[8, 4]
      betai9c4 <- beta[9, 4]
      betai10c4 <- beta[10, 4]
      bi1c1 <- coef(res)$X1[2]
      bi2c1 <- coef(res)$X2[2]
      bi3c1 <- coef(res)$X3[2]
      bi4c1 <- coef(res)$X4[2]
      bi5c1 <- coef(res)$X5[2]
      bi6c1 <- coef(res)$X6[2]
      bi7c1 <- coef(res)$X7[2]
      bi8c1 <- coef(res)$X8[2]
      bi9c1 <- coef(res)$X9[2]
      bi10c1 <- coef(res)$X10[2]
      bi1c2 <- coef(res)$X1[3]
      bi2c2 <- coef(res)$X2[3]
      bi3c2 <- coef(res)$X3[3]
      bi4c2 <- coef(res)$X4[3]
      bi5c2 <- coef(res)$X5[3]
      bi6c2 <- coef(res)$X6[3]
      bi7c2 <- coef(res)$X7[3]
      bi8c2 <- coef(res)$X8[3]
      bi9c2 <- coef(res)$X9[3]
      bi10c2 <- coef(res)$X10[3]
      bi1c3 <- coef(res)$X1[4]
      bi2c3 <- coef(res)$X2[4]
      bi3c3 <- coef(res)$X3[4]
      bi4c3 <- coef(res)$X4[4]
      bi5c3 <- coef(res)$X5[4]
      bi6c3 <- coef(res)$X6[4]
      bi7c3 <- coef(res)$X7[4]
      bi8c3 <- coef(res)$X8[4]
      bi9c3 <- coef(res)$X9[4]
      bi10c3 <- coef(res)$X10[4]
      bi1c4 <- coef(res)$X1[5]
      bi2c4 <- coef(res)$X2[5]
      bi3c4 <- coef(res)$X3[5]
      bi4c4 <- coef(res)$X4[5]
      bi5c4 <- coef(res)$X5[5]
      bi6c4 <- coef(res)$X6[5]
      bi7c4 <- coef(res)$X7[5]
      bi8c4 <- coef(res)$X8[5]
      bi9c4 <- coef(res)$X9[5]
      bi10c4 <- coef(res)$X10[5]
    }
    out <- tibble::tibble(numcat, alpha1, alpha2, alpha3, alpha4, alpha5,
                          alpha6, alpha7, alpha8, alpha9, alpha10,
                          a1, a2, a3, a4, a5,
                          a6, a7, a8, a9, a10,
                          betai1c1, betai2c1, betai3c1, betai4c1, betai5c1,
                          betai6c1, betai7c1, betai8c1, betai9c1, betai10c1,
                          betai1c2, betai2c2, betai3c2, betai4c2, betai5c2,
                          betai6c2, betai7c2, betai8c2, betai9c2, betai10c2,
                          betai1c3, betai2c3, betai3c3, betai4c3, betai5c3,
                          betai6c3, betai7c3, betai8c3, betai9c3, betai10c3,
                          betai1c4, betai2c4, betai3c4, betai4c4, betai5c4,
                          betai6c4, betai7c4, betai8c4, betai9c4, betai10c4,
                          bi1c1, bi2c1, bi3c1, bi4c1, bi5c1,
                          bi6c1, bi7c1, bi8c1, bi9c1, bi10c1,
                          bi1c2, bi2c2, bi3c2, bi4c2, bi5c2,
                          bi6c2, bi7c2, bi8c2, bi9c2, bi10c2,
                          bi1c3, bi2c3, bi3c3, bi4c3, bi5c3,
                          bi6c3, bi7c3, bi8c3, bi9c3, bi10c3,
                          bi1c4, bi2c4, bi3c4, bi4c4, bi5c4,
                          bi6c4, bi7c4, bi8c4, bi9c4, bi10c4)
    return(out)
  }
  #========================================================================
  # MAE and RMSE
  #========================================================================
  mymae <- function(data, dich = T, from = 1.6, to = 1.8, by = .001) {
    library(LaplacesDemon)
    d <- seq(from = from, to = to, by = by)
    mae <- bestd <-  vector("double", length(d))
    dat <- data
    for (i in 1:length(d)) {
      if (dich) {
        dat <- dat %>% mutate(ab1 = abs(alpha1 - a1 / d[i]),
                              ab2 = abs(alpha2 - a2 / d[i]),
                              ab3 = abs(alpha3 - a3 / d[i]),
                              ab4 = abs(alpha4 - a4 / d[i]),
                              ab5 = abs(alpha5 - a5 / d[i]),
                              ab6 = abs(alpha6 - a6 / d[i]),
                              ab7 = abs(alpha7 - a7 / d[i]),
                              ab8 = abs(alpha8 - a8 / d[i]),
                              ab9 = abs(alpha9 - a9 / d[i]),
                              ab10 = abs(alpha10 - a10 / d[i]),
                              ab1c1 = abs(betai1c1 + bi1c1 / d[i]),
                              ab2c1 = abs(betai2c1 + bi2c1 / d[i]),
                              ab3c1 = abs(betai3c1 + bi3c1 / d[i]),
                              ab4c1 = abs(betai4c1 + bi4c1 / d[i]),
                              ab5c1 = abs(betai5c1 + bi5c1 / d[i]),
                              ab6c1 = abs(betai6c1 + bi6c1 / d[i]),
                              ab7c1 = abs(betai7c1 + bi7c1 / d[i]),
                              ab8c1 = abs(betai8c1 + bi8c1 / d[i]),
                              ab9c1 = abs(betai9c1 + bi9c1 / d[i]),
                              ab10c1 = abs(betai10c1 + bi10c1 / d[i]))
        dat <- dat %>% mutate(absum = ab1 + ab2 + ab3 + ab4 + ab5 +
                                ab6 + ab7 + ab8 + ab9 + ab10 +
                                ab1c1 + ab2c1 + ab3c1 + ab4c1 + ab5c1 +
                                ab6c1 + ab7c1 + ab8c1 + ab9c1 + ab10c1)
        temp <- dat %>% dplyr::select(absum) %>% summarize_all(mean)
        mae[i] <- sqrt(temp / 20)
      } else {
        dat <- dat %>% mutate(ab1 = abs(alpha1 - a1 / d[i]),
                              ab2 = abs(alpha2 - a2 / d[i]),
                              ab3 = abs(alpha3 - a3 / d[i]),
                              ab4 = abs(alpha4 - a4 / d[i]),
                              ab5 = abs(alpha5 - a5 / d[i]),
                              ab6 = abs(alpha6 - a6 / d[i]),
                              ab7 = abs(alpha7 - a7 / d[i]),
                              ab8 = abs(alpha8 - a8 / d[i]),
                              ab9 = abs(alpha9 - a9 / d[i]),
                              ab10 = abs(alpha10 - a10 / d[i]),
                              ab1c1 = abs(betai1c1 + bi1c1 / d[i]),
                              ab2c1 = abs(betai2c1 + bi2c1 / d[i]),
                              ab3c1 = abs(betai3c1 + bi3c1 / d[i]),
                              ab4c1 = abs(betai4c1 + bi4c1 / d[i]),
                              ab5c1 = abs(betai5c1 + bi5c1 / d[i]),
                              ab6c1 = abs(betai6c1 + bi6c1 / d[i]),
                              ab7c1 = abs(betai7c1 + bi7c1 / d[i]),
                              ab8c1 = abs(betai8c1 + bi8c1 / d[i]),
                              ab9c1 = abs(betai9c1 + bi9c1 / d[i]),
                              ab10c1 = abs(betai10c1 + bi10c1 / d[i]),
                              ab1c2 = abs(betai1c2 + bi1c2 / d[i]),
                              ab2c2 = abs(betai2c2 + bi2c2 / d[i]),
                              ab3c2 = abs(betai3c2 + bi3c2 / d[i]),
                              ab4c2 = abs(betai4c2 + bi4c2 / d[i]),
                              ab5c2 = abs(betai5c2 + bi5c2 / d[i]),
                              ab6c2 = abs(betai6c2 + bi6c2 / d[i]),
                              ab7c2 = abs(betai7c2 + bi7c2 / d[i]),
                              ab8c2 = abs(betai8c2 + bi8c2 / d[i]),
                              ab9c2 = abs(betai9c2 + bi9c2 / d[i]),
                              ab10c2 = abs(betai10c2 + bi10c2 / d[i]),
                              ab1c3 = abs(betai1c3 + bi1c3 / d[i]),
                              ab2c3 = abs(betai2c3 + bi2c3 / d[i]),
                              ab3c3 = abs(betai3c3 + bi3c3 / d[i]),
                              ab4c3 = abs(betai4c3 + bi4c3 / d[i]),
                              ab5c3 = abs(betai5c3 + bi5c3 / d[i]),
                              ab6c3 = abs(betai6c3 + bi6c3 / d[i]),
                              ab7c3 = abs(betai7c3 + bi7c3 / d[i]),
                              ab8c3 = abs(betai8c3 + bi8c3 / d[i]),
                              ab9c3 = abs(betai9c3 + bi9c3 / d[i]),
                              ab10c3 = abs(betai10c3 + bi10c3 / d[i]),
                              ab1c4 = abs(betai1c4 + bi1c4 / d[i]),
                              ab2c4 = abs(betai2c4 + bi2c4 / d[i]),
                              ab3c4 = abs(betai3c4 + bi3c4 / d[i]),
                              ab4c4 = abs(betai4c4 + bi4c4 / d[i]),
                              ab5c4 = abs(betai5c4 + bi5c4 / d[i]),
                              ab6c4 = abs(betai6c4 + bi6c4 / d[i]),
                              ab7c4 = abs(betai7c4 + bi7c4 / d[i]),
                              ab8c4 = abs(betai8c4 + bi8c4 / d[i]),
                              ab9c4 = abs(betai9c4 + bi9c4 / d[i]),
                              ab10c4 = abs(betai10c4 + bi10c4 / d[i]))
        dat <- dat %>% mutate(absum = ab1 + ab2 + ab3 + ab4 + ab5 +
                                ab6 + ab7 + ab8 + ab9 + ab10 +
                                ab1c1 + ab2c1 + ab3c1 + ab4c1 + ab5c1 +
                                ab6c1 + ab7c1 + ab8c1 + ab9c1 + ab10c1 +
                                ab1c2 + ab2c2 + ab3c2 + ab4c2 + ab5c2 +
                                ab6c2 + ab7c2 + ab8c2 + ab9c2 + ab10c2 +
                                ab1c3 + ab2c3 + ab3c3 + ab4c3 + ab5c3 +
                                ab6c3 + ab7c3 + ab8c3 + ab9c3 + ab10c3 +
                                ab1c4 + ab2c4 + ab3c4 + ab4c4 + ab5c4 +
                                ab6c4 + ab7c4 + ab8c4 + ab9c4 + ab10c4)
        temp <- dat %>% dplyr::select(absum) %>% summarize_all(mean)
        mae[i] <- sqrt(temp / 50)
      }
    }
    mae <- unlist(mae)
    bestd <- d[which(mae == min(mae))]
    out <- list(mae, bestd)
    return(out)
  }
  myrmse <- function(data, dich = T, from = 1.6, to = 1.8, by = .001) {
    library(LaplacesDemon)
    d <- seq(from = from, to = to, by = by)
    rmse <- bestd <-  vector("double", length(d))
    dat <- data
    for (i in 1:length(d)) {
      if (dich) {
        dat <- dat %>% mutate(sq1 = (alpha1 - a1 / d[i]) ^ 2,
                              sq2 = (alpha2 - a2 / d[i]) ^ 2,
                              sq3 = (alpha3 - a3 / d[i]) ^ 2,
                              sq4 = (alpha4 - a4 / d[i]) ^ 2,
                              sq5 = (alpha5 - a5 / d[i]) ^ 2,
                              sq6 = (alpha6 - a6 / d[i]) ^ 2,
                              sq7 = (alpha7 - a7 / d[i]) ^ 2,
                              sq8 = (alpha8 - a8 / d[i]) ^ 2,
                              sq9 = (alpha9 - a9 / d[i]) ^ 2,
                              sq10 = (alpha10 - a10 / d[i]) ^ 2,
                              sq1c1 = (betai1c1 + bi1c1 / d[i]) ^2,
                              sq2c1 = (betai2c1 + bi2c1 / d[i]) ^2,
                              sq3c1 = (betai3c1 + bi3c1 / d[i]) ^2,
                              sq4c1 = (betai4c1 + bi4c1 / d[i]) ^2,
                              sq5c1 = (betai5c1 + bi5c1 / d[i]) ^2,
                              sq6c1 = (betai6c1 + bi6c1 / d[i]) ^2,
                              sq7c1 = (betai7c1 + bi7c1 / d[i]) ^2,
                              sq8c1 = (betai8c1 + bi8c1 / d[i]) ^2,
                              sq9c1 = (betai9c1 + bi9c1 / d[i]) ^2,
                              sq10c1 = (betai10c1 + bi10c1 / d[i]) ^2)
        dat <- dat %>% mutate(sumsqre = sq1 + sq2 + sq3 + sq4 + sq5 +
                                sq6 + sq7 + sq8 + sq9 + sq10 +
                                sq1c1 + sq2c1 + sq3c1 + sq4c1 + sq5c1 +
                                sq6c1 + sq7c1 + sq8c1 + sq9c1 + sq10c1)
        temp <- dat %>% dplyr::select(sumsqre) %>% summarize_all(mean)
        rmse[i] <- sqrt(temp / 20)
      } else {
        dat <- dat %>% mutate(sq1 = (alpha1 - a1 / d[i]) ^ 2,
                              sq2 = (alpha2 - a2 / d[i]) ^ 2,
                              sq3 = (alpha3 - a3 / d[i]) ^ 2,
                              sq4 = (alpha4 - a4 / d[i]) ^ 2,
                              sq5 = (alpha5 - a5 / d[i]) ^ 2,
                              sq6 = (alpha6 - a6 / d[i]) ^ 2,
                              sq7 = (alpha7 - a7 / d[i]) ^ 2,
                              sq8 = (alpha8 - a8 / d[i]) ^ 2,
                              sq9 = (alpha9 - a9 / d[i]) ^ 2,
                              sq10 = (alpha10 - a10 / d[i]) ^ 2,
                              sq1c1 = (betai1c1 + bi1c1 / d[i]) ^2,
                              sq2c1 = (betai2c1 + bi2c1 / d[i]) ^2,
                              sq3c1 = (betai3c1 + bi3c1 / d[i]) ^2,
                              sq4c1 = (betai4c1 + bi4c1 / d[i]) ^2,
                              sq5c1 = (betai5c1 + bi5c1 / d[i]) ^2,
                              sq6c1 = (betai6c1 + bi6c1 / d[i]) ^2,
                              sq7c1 = (betai7c1 + bi7c1 / d[i]) ^2,
                              sq8c1 = (betai8c1 + bi8c1 / d[i]) ^2,
                              sq9c1 = (betai9c1 + bi9c1 / d[i]) ^2,
                              sq10c1 = (betai10c1 + bi10c1 / d[i]) ^2,
                              sq1c2 = (betai1c2 + bi1c2 / d[i]) ^2,
                              sq2c2 = (betai2c2 + bi2c2 / d[i]) ^2,
                              sq3c2 = (betai3c2 + bi3c2 / d[i]) ^2,
                              sq4c2 = (betai4c2 + bi4c2 / d[i]) ^2,
                              sq5c2 = (betai5c2 + bi5c2 / d[i]) ^2,
                              sq6c2 = (betai6c2 + bi6c2 / d[i]) ^2,
                              sq7c2 = (betai7c2 + bi7c2 / d[i]) ^2,
                              sq8c2 = (betai8c2 + bi8c2 / d[i]) ^2,
                              sq9c2 = (betai9c2 + bi9c2 / d[i]) ^2,
                              sq10c2 = (betai10c2 + bi10c2 / d[i]) ^2,
                              sq1c3 = (betai1c3 + bi1c3 / d[i]) ^2,
                              sq2c3 = (betai2c3 + bi2c3 / d[i]) ^2,
                              sq3c3 = (betai3c3 + bi3c3 / d[i]) ^2,
                              sq4c3 = (betai4c3 + bi4c3 / d[i]) ^2,
                              sq5c3 = (betai5c3 + bi5c3 / d[i]) ^2,
                              sq6c3 = (betai6c3 + bi6c3 / d[i]) ^2,
                              sq7c3 = (betai7c3 + bi7c3 / d[i]) ^2,
                              sq8c3 = (betai8c3 + bi8c3 / d[i]) ^2,
                              sq9c3 = (betai9c3 + bi9c3 / d[i]) ^2,
                              sq10c3 = (betai10c3 + bi10c3 / d[i]) ^2,
                              sq1c4 = (betai1c4 + bi1c4 / d[i]) ^2,
                              sq2c4 = (betai2c4 + bi2c4 / d[i]) ^2,
                              sq3c4 = (betai3c4 + bi3c4 / d[i]) ^2,
                              sq4c4 = (betai4c4 + bi4c4 / d[i]) ^2,
                              sq5c4 = (betai5c4 + bi5c4 / d[i]) ^2,
                              sq6c4 = (betai6c4 + bi6c4 / d[i]) ^2,
                              sq7c4 = (betai7c4 + bi7c4 / d[i]) ^2,
                              sq8c4 = (betai8c4 + bi8c4 / d[i]) ^2,
                              sq9c4 = (betai9c4 + bi9c4 / d[i]) ^2,
                              sq10c4 = (betai10c4 + bi10c4 / d[i]) ^2)
        dat <- dat %>% mutate(sumsqre = sq1 + sq2 + sq3 + sq4 + sq5 +
                                sq6 + sq7 + sq8 + sq9 + sq10 +
                                sq1c1 + sq2c1 + sq3c1 + sq4c1 + sq5c1 +
                                sq6c1 + sq7c1 + sq8c1 + sq9c1 + sq10c1 +
                                sq1c2 + sq2c2 + sq3c2 + sq4c2 + sq5c2 +
                                sq6c2 + sq7c2 + sq8c2 + sq9c2 + sq10c2 +
                                sq1c3 + sq2c3 + sq3c3 + sq4c3 + sq5c3 +
                                sq6c3 + sq7c3 + sq8c3 + sq9c3 + sq10c3 +
                                sq1c4 + sq2c4 + sq3c4 + sq4c4 + sq5c4 +
                                sq6c4 + sq7c4 + sq8c4 + sq9c4 + sq10c4)
        temp <- dat %>% dplyr::select(sumsqre) %>% summarize_all(mean)
        rmse[i] <- sqrt(temp / 50)
      }
    }
    rmse <- unlist(rmse)
    bestd <- d[which(rmse == min(rmse))]
    out <- list(rmse, bestd)
    return(out)
  }
  #========================================================================
  # Loop
  #========================================================================
  for (rep_set in rep_sets) {
    tictoc::tic()
    print(paste("Starting number of item ", numcat, ", rep_set", rep_set))
    filename <- paste0(name, numcat, "-", rep_set, ".csv")
    if (!file.exists(filename)) {
      out <- matrix(vector("double", length(rep_per_sets) * out_ncol), ncol = out_ncol)
      colnames(out) <- c("numcat", "rep_set", "rep",
                         paste0("Model", 1:10, stats[1]),
                         paste0("Model", 1:10, stats[2]))
      for (rep in rep_per_sets) {
        cat(name, ": ", numcat, ": rep set: ", rep_set, "rep: ", rep)
        out[rep, 1] <- numcat
        out[rep, 2] <- rep_set
        out[rep, 3] <- rep
        for (model in models) {
          cat(name, ": ", numcat, ": rep set: ", rep_set, "rep: ", rep, "model: ", model)
          genout <- generate(numcat, model, rep_set, rep)
          analout <- analyze(numcat, genout)
          mae <- mymae(analout, dich, optmae, optmae)[[1]]
          rmse <- myrmse(analout, dich, optrmse, optrmse)[[1]]
          out[rep, 3 + model] <- mae
          out[rep, 13 + model] <- rmse
        }
      } # end of for (rep in rep_per_set)
      readr::write_csv(data.frame(out), file = filename)
      print(out)
      tictoc::toc()
    } # end of if (!file.exists(filename))
  } # end of  for (rep_set in rep_sets)
  #========================================================================
  # Summarize
  #========================================================================
  dich_data <- 1:20 %>%
    purrr::map(~ read_csv(file = paste0(name, 2, "-", .x, ".csv"))) %>%
    dplyr::bind_rows()
  poly_data <- 1:20 %>%
    purrr::map(~ read_csv(file = paste0(name, 5, "-", .x, ".csv"))) %>%
    dplyr::bind_rows()
  out <- dplyr::bind_rows(dich_data, poly_data) %>%
    group_by(numcat) %>%
    summarize_all(mean) %>%
    select(-rep_set, -rep)
  write_csv(out, "study6.csv")
  return(out)
} # end of function
