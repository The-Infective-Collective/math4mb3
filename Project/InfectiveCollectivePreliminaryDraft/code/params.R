initfun <- function(param) {
    with(param,{
        epsilon <- mu/(mu+gamma)
        
        list(
            S=pop/R0,
            I=epsilon*(1-1/R0)*pop,
            R=pop-(1-R0-epsilon*(1-1/R0))
        )
    })
}

base.params <- list(
    R0=17,
    pop=1e6,
    b1=0.08,
    gamma=365/13,
    mu=0.02,
    dt=1/365,
    nsteps=365*100
)

base.init <- initfun(base.params)

base.M <- matrix(1/2, 2, 2)
