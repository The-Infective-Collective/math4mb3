params <- list()
params <- within(params, {
  
  # RNG seed
  seed <- 0
  
  # time step length (1 day)
  dt <- 1/365
  
  # number of years
  nyears <- 100
  
  # number of time steps
  nsteps <- nyears/dt
  
  # number of spatial patches (2 for most purposes, more for exploring mixing structure)
  mpatches <- 2
  
  # birth/death rate
  mu <- 1/70
  
  # recovery rate
  gamma <- 12/365 
  
  # R0
  R0 <- 17
  
  # population size in each patch
  pop <- 1e6
  
  # mean transmission rate <b_0>
  b_0 <- 0.6
  
  # seasonal amplitude
  b_1 <- 0.1
  
  # 2*pi divided by period length??
  a <- 2.0
  
  # connectivity matrix
  conmat <- matrix(rep(mpatches^(-1), mpatches^2), nrow = mpatches, ncol = mpatches)
  
  # initial conditions, 
  init <- within(list(), {
    S <- round(pop/R0)
    I <- round(mu/(mu + gamma)*(1 - 1/R0)*pop)
    R <- pop-S-I
  })
})

