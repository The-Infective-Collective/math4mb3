library(Rcpp)
library(tidyverse,stringr)
sourceCpp("SIRmodel_mix.cpp")
source("params.R")

# change number of patches
params$mpatches <- 6

# set seed
set.seed(params$seed)

# function to run model and return data frame
run_SIR <- function(p) {
  # Run the model
  raw.result <- SIRmodel(p)
  
  # Initialize result data frame
  result <- data_frame(t = raw.result$time)
  
  # Pull results
  S <- raw.result$S
  I <- raw.result$I
  R <- raw.result$R
  
  # Change column names
  colnames(S) <- str_c("S", 1:p$mpatches)
  colnames(I) <- str_c("I", 1:p$mpatches)
  colnames(R) <- str_c("R", 1:p$mpatches)
  
  # Combine data frames 
  df <- cbind(result, S, I, R) 
  
  # melt for easy plotting
  melt(df, c("t"), str_c("I", 1:params$mpatches), variable.name = "compartment", value.name = "prevalence")
}

# set plotting theme for ggplot
theme_set(theme_bw())

# run model
result <- run_SIR(params)
result

# plot
ggplot(result, aes(t, prevalence, colour = compartment)) +
  geom_line()

