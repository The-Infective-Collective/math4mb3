# convenient functions

# function to run model and return data frame
run_SIR <- function(p, init, m, seas) {
  # get number of patches
  n <- nrow(m)
  
  # Run the model
  raw.result <- SIRmodel_npatch(p, init, m, seas)
  
  # Initialize result data frame
  result <- data_frame(t = raw.result$time)
  
  # Pull results
  S <- raw.result$S
  I <- raw.result$I
  R <- raw.result$R
  
  # Change column names
  colnames(S) <- str_c("S", 1:n)
  colnames(I) <- str_c("I", 1:n)
  colnames(R) <- str_c("R", 1:n)
  
  # Combine data frames 
  cbind(result, S, I, R) 
}

# coherence metric calculation
coherence_calc <- function(x) {
  norm(x - rep(mean(x), length(x)), type = "2")
}

# coherence metric scaled to the size of the epidemic
coherence_calc_scaled <- function(x) {
  norm(x - rep(mean(x), length(x)), type = "2") / mean(x)
}

# Generate initial conditions
## generate by adding noise to the endemic equilibrium
r.init <- function(m) { 
  initial <- list(
    S = round(rep(base.init$S, m) + runif(m, -base.init$S, base.init$S)),
    I = round(rep(base.init$I, m) + runif(m, - base.init$I, 30*base.init$I))
  )
  initial$R <- rep(base.params$pop, m) - initial$S - initial$I
  initial
}

## generate by sampling uniform 
r.init.unif <- function(m) { 
  initial <- list(S = runif(m, 0, base.params$pop))
  initial$I <- runif(m, 0, base.params$pop - initial$S)
  initial$R <- rep(base.params$pop, m) - initial$S - initial$I
  initial
}

r.init.science <- function(m) {
    initial <- list(S = runif(m, 0, 0.1*base.params$pop))
    initial$I <- runif(m, 0, 0.0001*base.params$pop)
    initial$R <- rep(base.params$pop, m) - initial$S - initial$I
    initial
}
