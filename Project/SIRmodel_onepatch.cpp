#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List SIRmodel(List params, List init) {
    
    // Number of time steps
    int nsteps = as<int>(params["nsteps"]);
    
    // Initalize matrices to hold results (time steps as rows, columns as patches) - makes matrix of zeroes
    NumericVector SS(nsteps);
    NumericVector II(nsteps);
    NumericVector RR(nsteps);
    
    // Set initial conditions
    SS[0] = init["S"];
    II[0] = init["I"];
    RR[0] = init["R"];
    
    
    // Initialize time vector
    NumericVector time(nsteps);
    
    time[0] = 0;
    
    NumericVector incidence(nsteps);
    
    // Pull out parameter values
    double dt = params["dt"]; 
    double mu = params["mu"];
    double gamma = params["gamma"];
    double recoveryrate = 1 - exp(-(gamma+mu)*dt);
    double deathrate = 1 - exp(-mu*dt);
    
    double b0 = params["b0"];
    double b1 = params["b1"];
    double a = params["a"];
    
    double beta;
    double leaveS, leaveI, leaveR;
    double enterS, enterI, enterR;
    
    // Iterate over time steps
    for (int istep = 0; istep < (nsteps-1); istep++) {
        beta = b0 * (1 + b1 * cos (a * time[istep]));
        
        // The proportion that leave the susceptible class
        leaveS = (1-exp(-(beta*II[istep]+mu)*dt))*SS[istep];
        leaveI = recoveryrate*II[istep];
        leaveR = deathrate*RR[istep];
        
        enterI = beta*II[istep]/(beta*II[istep]+mu) * leaveS;
        enterR = gamma/(gamma+mu) * leaveI;
        enterS = leaveS - enterI + leaveI - enterR + leaveR;
        
        // Update state for next time step according to model
        SS[istep+1] = SS(istep) + enterS - leaveS;  
        II[istep+1] = II(istep) + enterI - leaveI;
        RR[istep+1] = RR(istep) + enterR - leaveR;
        
        incidence[istep] = enterI;
        
        // Time in fraction of year
        time[istep+1] = (istep+1)*dt;
    }
    
    // Initialize a list to return results
    List ret;
    
    // Put results in a list
    ret["time"] = time;
    ret["S"] = SS;
    ret["I"] = II;
    ret["R"] = RR;
    ret["incidence"] = incidence;
    
    return ret;
};
