library("deSolve")

## Vector Field for SIR model
SIR.vector.field <- function(t, vars, parms=c(beta=betavals,gamma=0.25)) {
  with(as.list(c(parms, vars)), {
    dx <- -beta*x*y # dS/dt
    dy <- beta*x*y - gamma*y # dI/dt
    vec.fld <- c(dx=dx, dy=dy)
    return(list(vec.fld)) # ode() requires a list
  })
}

## install tikzDevice
#library("tikzDevice")
#tikz("SIRsolns_graph.tex",standAlone=TRUE)

## Draw solution

## Let S0 rep initial proportion of susceptible individuals
## Let I0 rep initial proportion of infected individuals
## Let beta rep transmission rate
## Let gamma rep recovery rate
draw.soln <- function(ic=c(x=1,y=0), tmax=1,
                      times=seq(0,tmax,by=tmax/500),
                      func, parms, ... ) {
  soln <- ode(ic, times, func, parms)
  lines(times, soln[,"y"], col=i, lwd=3, ... )
}

## Plot solutions of the SIR model
tmax <- 150 # end time for numerical integration of the ODE

## draw box for plot:
plot(0,0,xlim=c(0,tmax),ylim=c(0,0.5),
     type="n",xlab="Time (days)",ylab="Prevalence",las=1)

## initial conditions:
I0 <- 0.001
S0 <- 1 - I0
R0 <- 0

## parameters
R_0vals <- c(1.2,1.5,1.8,2,3,4)
betavals <- 0.25*R_0vals
gamma <- 0.25

## draw solutions for several values of parameter beta:
R_0vals <- c(1.2,1.5,1.8,2,3,4)
betavals <- 0.25*R_0vals
for (i in 1:length(betavals)) {
  draw.soln(ic=c(x=S0,y=I0), tmax=tmax,
            func=SIR.vector.field,
            parms=c(beta=betavals[i],gamma=gamma),
            lty=i
  )
}
legend(115, 0.5, title="Reproductive number",
       legend=c("1.2", "1.5", "1.8", "2", "3", "4"),
       col=c("black", "red", "green", "blue", "cyan","magenta"),
       lty=c(1:6), lwd=2, cex=0.85, seg.len=4)