study5 <- function(start_con = 1, end_con = 10,
                start_mod = 1, end_mod = 13,
                start_D = 1.65, end_D = 1.75,
                start_rep = 1, end_rep = 10) {
  library(tidyverse)
  library(LaplacesDemon)
  library(semTools)
  library(psych)
  library(tictoc)
  numcat <- c(2, 5)
  dist <- c(1, 2, 3, 4, 5) # 1 normal, 2 skewed, 3 platykurtic, 4 leptokurtic, 5 severe
  conditions <- tidyr::crossing(numcat, dist)
  colnames(conditions) <- c("numcat", "dist")
  condition_numbers <- start_con:end_con
  models <- start_mod:end_mod
  D <- seq(start_D, end_D, length.out = 11)
  reps <- start_rep:end_rep
  name <- "study5-"
  #========================================================================
  # Analyze
  #========================================================================
  analyze <- function(conditions, condition_number, rep, what = "L", D = 1.702) {
    set.seed(10 * condition_number + rep)
    numcat <- as.integer(conditions[condition_number, 1])
    dist <- as.integer(conditions[condition_number, 2])
    n <- 10^7
    k <- 6
    alpha <- c(.5, .7, .9, 1.1, 1.3, 1.5)
    if (numcat == 2) {
      beta <- 0
    } else {
      beta <- c(-.9, -.3, .3, .9)
    }
    skew <- 0
    kurt <- 0
    if (dist == 2) {
      skew <- 2
    } else if (dist == 3) {
      kurt <-  -1.2
    } else if (dist == 4) {
      kurt <- 7
    } else if (dist == 5) {
      skew = 3
      kurt = 21
    }
    theta <- mvrnonnorm(n, c(0,0), matrix(c(1, 0, 0, 1), 2, 2),
                        skewness = c(skew, 0),
                        kurtosis = c(kurt, 0))[, 1]
    # Factor Analysis model
    if (what == "F") {
      load <- alpha / sqrt(1 + alpha ^ 2)
      if (numcat == 2) {
        thresh <- vector("double", k)
        for (j in 1:k) {
          thresh[j] <- beta / sqrt(1 + alpha[j] ^ 2)
        }
      } else {
        thresh <- matrix(vector("double", k * (numcat - 1)), nrow = k)
        for (j in 1:k) {
          for (l in 1:(numcat - 1)) {
            thresh[j, l] <- beta[l] / sqrt(1 + alpha[j] ^ 2)
          }
        }
      }
      e <- matrix(rnorm(n * k), nrow = n, ncol = k)
      score_con <- theta %*% t(load) + e %*% sqrt(diag(k) - diag(diag(load %*% t(load))))
      parallel_e <- matrix(rnorm(n * k), nrow = n, ncol = k)
      parallel_con <-theta %*% t(load) + parallel_e %*% sqrt(diag(k) - diag(diag(load %*% t(load))))
      score_ord <- matrix(vector("double", n * k), nrow = n)
      parallel_ord <- matrix(vector("double", n * k), nrow = n)
      for (i in 1:n) {
        for (j in 1:k) {
          if (numcat == 2) {
            score_ord[i, j] <- ifelse(score_con[i, j] < thresh[j], 0, 1)
            parallel_ord[i, j] <- ifelse(parallel_con[i, j] < thresh[j], 0, 1)
          } else {
            score_ord[i, j] <- ifelse(score_con[i, j] < thresh[j, 1], 0,
                                      ifelse(score_con[i, j] < thresh[j, 2], 1,
                                             ifelse(score_con[i, j] < thresh[j, 3], 2,
                                                    ifelse(score_con[i, j] < thresh[j, 4], 3, 4))))
            parallel_ord[i, j] <- ifelse(parallel_con[i, j] < thresh[j, 1], 0,
                                         ifelse(parallel_con[i, j] < thresh[j, 2], 1,
                                                ifelse(parallel_con[i, j] < thresh[j, 3], 2,
                                                       ifelse(parallel_con[i, j] < thresh[j, 4], 3, 4))))
          }
        }
      }
      out <- sum(cov(score_ord, parallel_ord)) / sqrt(sum(var(score_ord)) * sum(var(parallel_ord)))
    } else {
      # IRT model
      true_score <- vector("double", n)
      item_score <- matrix(vector("double", n * k), nrow = n)
      if (numcat == 2) {
        among <- c(0, 1)
        p <- matrix(vector("double", n * k), nrow = n)
        for (i in 1:n) {
          for (j in 1:k) {
            if (what == "N") {
              p[i, j] <- pnorm(alpha[j] * theta[i] - beta)
            } else {
              p[i, j] <- invlogit(D * alpha[j] * theta[i] - D * beta)
            }
            true_score[i] <- true_score[i] + p[i, j]
            item_score[i, j] <- sample(among, size = 1, replace = TRUE,
                                       prob = c(1 - p[i, j],
                                                p[i, j]))
          }
        }
      } else {
        among <- c(0, 1, 2, 3, 4)
        p <- array(vector("double", n * k * (numcat - 1)),
                   dim = c(n, k, (numcat - 1)))
        for (i in 1:n) {
          for (j in 1:k) {
            for (l in 1:(numcat - 1)) {
              if (what == "N") {
                p[i, j, l] <- pnorm(alpha[j] * theta[i] - beta[l])
              } else {
                p[i, j, l] <- invlogit(D * alpha[j] * theta[i] - D * beta[l])
              }
              true_score[i] <- true_score[i] + p[i, j, l]
            }
            item_score[i, j] <- sample(among, size = 1, replace = TRUE,
                                       prob = c(1 - p[i, j, 1],
                                                p[i, j, 1] - p[i, j, 2],
                                                p[i, j, 2] - p[i, j, 3],
                                                p[i, j, 3] - p[i, j, 4],
                                                p[i, j, 4]))
          }
        } # for (i in 1:n) {
      } # else: if (numcat == 2) {
      out <- sum(var(true_score)) / sum(var(item_score))
    } # } else: if (what == "F") {
    return(out)
  }
  #========================================================================
  # Loop
  #========================================================================
  for (condition_number in condition_numbers) {
    condition <- conditions[condition_number, ]
    print(condition)
    for (model in models) {
      if (model == 1 | model == 2) {
        Dname <- "NA"
      } else {
        Dname <- paste0(D[model - 2])
      }
      for (rep in reps) {
        print(paste0("Starting model number", model, " rep: ", rep))
        tictoc::tic()
        filename <- paste0(name, condition_number, "M", model, "r", rep, "-", Dname, ".csv")
        if (!file.exists(filename)) {
          if (model == 1) {
            out <- analyze(conditions, condition_number, rep,  "F")
          } else if (model == 2) {
            out <- analyze(conditions, condition_number, rep, "N")
          } else {
            out <- analyze(conditions, condition_number, rep, "L", D[model - 2])
          }
          readr::write_csv(data.frame(out), file = filename)
          tictoc::toc()
      }
      } # end of if (!file.exists(filename))
    } # end of  for (rep_set in rep_sets)
  } # end of for (condition_number in condition_numbers)
  #========================================================================
  # Summarize
  #========================================================================
  cons <- 10
  models <- 13
  reps <- 10
  D <- seq(from = 1.65, to = 1.75, by = .01)
  read_array <- array(vector("double", cons * models * reps), dim = c(cons, models, reps))
  read_df <- data.frame(matrix(vector("double", cons * models), nrow = cons))
  for (i in 1:cons) {
    for (j in 1:models) {
      tempsum <- 0
      for (l in 1:reps) {
        filename<- paste0("study5-", i, "M", j, "r", l, "-")
        if (j == 1 | j == 2) {
          filename<- paste0(filename, "NA")
        } else {
          filename<- paste0(filename, D[j - 2])
        }
        filename<- paste0(filename, ".csv")
        read_array[i, j, l] <- as.double(read.csv(filename))
        tempsum <- tempsum + read_array[i, j, l]
      }
      read_df[i, j] <- tempsum / reps
    }
  }
  colnames(read_df) <- c("FA", "NO", "D165", "D166", "D167", "D168", "D169",
                         "D170", "D171", "D172", "D173", "D174", "D175")
  rownames(read_df) <- c("normal2", "skewed2", "platy2", "lepto2", "severe2", 
                         "normal5", "skewed5", "platy5", "lepto5", "severe5")
  read_df <- read_df %>% select(-D165, -D174, -D175) %>% 
    relocate(NO, .after = D173) %>% 
    mutate_all(~ .x - NO) %>% 
    select(-NO)  
  write_csv(read_df, "study5.csv")
  return(read_df)
}
