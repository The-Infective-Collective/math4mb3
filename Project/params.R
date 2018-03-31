initfun <- function(param,
                    n.patch=1,
                    round=FALSE) {
    with(param,{
        epsilon <- mu/(mu+gamma)
        
        ll <- list(
            S=rep(pop/R0, n.patch),
            I=rep(epsilon*(1-1/R0)*pop, n.patch)
        )
        
        if (round) ll <- lapply(ll, round)
        
        ll$R <- pop - ll$S - ll$I
        
        ll
    })
}

base.params <- list(
    R0=17,
    pop=1e6,
    b1=0.25,
    gamma=365/13,
    mu=0.02,
    dt=1/365,
    nsteps=365*100
)

base.init <- initfun(base.params)

base.M <- matrix(1/2, 2, 2)
