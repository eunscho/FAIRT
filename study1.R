study1 <- function() {
  #========================================================================
  # minimizing the maximum error
  #========================================================================
  d <- seq(from = 1.6, to = 1.8, by = .0001)
  errors <- vector("double", length(d))
  # ogive function
  for (i in 1:length(d)) {
    fun1 <- function(x) pnorm(x) - plogis(x, scale = 1/d[i])
    fun2 <- function(x) -pnorm(x) + plogis(x, scale = 1/d[i])
    tail <- abs(optimize(fun1, interval = c(-4, 0))$objective)
    shoulder <- abs(optimize(fun2, interval = c(-4, 0))$objective)
    errors[i] <- max(tail, shoulder)
  }
  ogive_minimax <- d[which(errors == min(errors))]
  # pdf
  for (i in 1:length(d)) {
    fun1 <- function(x) dnorm(x) - dlogis(x, scale = 1/d[i])
    fun2 <- function(x) -dnorm(x) + dlogis(x, scale = 1/d[i])
    tail <- abs(optimize(fun1, interval = c(-4, -.5))$objective)
    shoulder <- abs(optimize(fun2, interval = c(-4, 0))$objective)
    center <- dlogis(0, scale = 1/d[i]) - dnorm(0)
    errors[i] <- max(tail, shoulder, center)
  }
  pdf_minimax <- d[which(errors == min(errors))]
  #========================================================================
  # minimizing the sum of errors
  #========================================================================
  use <- function(ogive = TRUE, squared = TRUE) {
    library(LaplacesDemon)
    for (i in 1:length(d)) {
      if (squared & ogive) {
        func <- function(x) {(invlogit(d[i] * x) - pnorm(x))^2}  
      } else if (!squared & ogive) {
        func <- function(x) {abs(invlogit(d[i] * x) - pnorm(x))}   
      } else if (squared & !ogive) {
        func <- function(x) {(dlogis(x, scale = 1 / d[i]) - dnorm(x))^2}  
      } else {
        func <- function(x) {abs(dlogis(x, scale = 1 / d[i]) - dnorm(x))}   
      }
      errors[i] <- integrate(func, -Inf, Inf)
    }
    errors <- unlist(errors)
    return(d[which(errors == min(errors))])
  }
  ogive_squared <- use(T, T)
  ogive_absolute <- use(T, F)
  pdf_squared <- use(F, T)
  pdf_absolute <- use(F, F)
  #========================================================================
  # reporting output
  #========================================================================
  out <- data.frame(ogive_minimax, ogive_squared, ogive_absolute, 
                    pdf_minimax, pdf_squared, pdf_absolute)
  return(out)
}

