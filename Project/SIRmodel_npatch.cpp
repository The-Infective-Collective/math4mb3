#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List SIRmodel_npatch(List params, List init, NumericMatrix m,
                     Rcpp::Function betafun) {
    
    // Number of time steps
    int nsteps = as<int>(params["nsteps"]);
    
    // Number of patches
    int mpatches = m.ncol();
    
    // Initalize matrices to hold results (time steps as rows, columns as patches) - makes matrix of zeroes
    NumericMatrix SS(nsteps, mpatches);
    NumericMatrix II(nsteps, mpatches);
    NumericMatrix RR(nsteps, mpatches);
    
    // Get the initial state values
    NumericVector Sinit = init["S"];
    NumericVector Iinit = init["I"];
    NumericVector Rinit = init["R"];
    
    // Set initial conditions in the first row
    for (int i = 0; i < mpatches; i++) {
        SS(0, i) = Sinit[i];
        II(0, i) = Iinit[i];
        RR(0, i) = Rinit[i];
    }
    
    // Initialize time vector
    NumericVector time(nsteps);
    
    // Pull out parameter values
    double dt = params["dt"]; 
    double mu = params["mu"];
    double gamma = params["gamma"];
    double recoveryrate = 1 - exp(-(gamma+mu)*dt);
    double deathrate = 1 - exp(-mu*dt);
    
    double beta;
    double leaveS, leaveI, leaveR;
    double enterS, enterI, enterR;
    
    // Iterate over time steps
    for (int istep = 0; istep < (nsteps-1); istep++) {
        
        beta = Rcpp::as<double>(betafun(time[istep], params));
        
        // Iterate over patches
        for (int jpatch = 0; jpatch < mpatches; jpatch++) {
            
            // Get current S, I, and R states
            double iS1 = SS(istep, jpatch);  
            double iI1 = II(istep, jpatch);
            double iR1 = RR(istep, jpatch);
            
            // Initialize variable to hold the sum of (connectivity)*(# infected)
            double transsum = 0;
            
            // Iterate over all patches
            for (int i=0; i < (mpatches); i++) {
                // Calculate sum of (connectivity)*(# infected)
                transsum += m(jpatch, i)*II(istep,i);
            }
            
            leaveS = (1-exp(-(beta*transsum+mu)*dt))*iS1;
            leaveI = recoveryrate*iI1;
            leaveR = deathrate*iR1;
            
            enterI = beta*transsum/(beta*transsum+mu) * leaveS;
            enterR = gamma/(gamma+mu) * leaveI;
            enterS = leaveS - enterI + leaveI - enterR + leaveR;
            
            // Update state for next time step according to model
            SS(istep+1, jpatch) = iS1 + enterS - leaveS;  
            II(istep+1, jpatch) = iI1 + enterI - leaveI;
            RR(istep+1, jpatch) = iR1 + enterR - leaveR;
        }
        
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
    
    return ret;
};

// [[Rcpp::export]]
List SIRmodel_npatch_stochastic(List params, List init, NumericMatrix m,
                                Rcpp::Function betafun) {
    
    // Number of time steps
    int nsteps = as<int>(params["nsteps"]);
    
    // Number of patches
    int mpatches = m.ncol();
    
    // Initalize matrices to hold results (time steps as rows, columns as patches) - makes matrix of zeroes
    NumericMatrix SS(nsteps, mpatches);
    NumericMatrix II(nsteps, mpatches);
    NumericMatrix RR(nsteps, mpatches);
    
    // Get the initial state values
    NumericVector Sinit = init["S"];
    NumericVector Iinit = init["I"];
    NumericVector Rinit = init["R"];
    
    // Set initial conditions in the first row
    for (int i = 0; i < mpatches; i++) {
        SS(0, i) = Sinit[i];
        II(0, i) = Iinit[i];
        RR(0, i) = Rinit[i];
    }
    
    // Initialize time vector
    NumericVector time(nsteps);
    
    // Pull out parameter values
    double dt = params["dt"]; 
    double mu = params["mu"];
    double gamma = params["gamma"];
    double recoveryrate = 1 - exp(-(gamma+mu)*dt);
    double deathrate = 1 - exp(-mu*dt);
    
    double beta;
    double leaveS, leaveI, leaveR;
    double enterS, enterI, enterR;
    
    // Iterate over time steps
    for (int istep = 0; istep < (nsteps-1); istep++) {
        
        beta = Rcpp::as<double>(betafun(time[istep], params));
        
        // Iterate over patches
        for (int jpatch = 0; jpatch < mpatches; jpatch++) {
            
            // Get current S, I, and R states
            double iS1 = SS(istep, jpatch);  
            double iI1 = II(istep, jpatch);
            double iR1 = RR(istep, jpatch);
            
            // Initialize variable to hold the sum of (connectivity)*(# infected)
            double transsum = 0;
            
            // Iterate over all patches
            for (int i=0; i < (mpatches); i++) {
                // Calculate sum of (connectivity)*(# infected)
                transsum += m(jpatch, i)*II(istep,i);
            }
            
            leaveS = R::rbinom(iS1, (1-exp(-(beta*transsum+mu)*dt)));
            leaveI = R::rbinom(iI1, recoveryrate);
            leaveR = R::rbinom(iR1, deathrate);
            
            enterI = R::rbinom(leaveS, beta*transsum/(beta*transsum+mu));
            enterR = R::rbinom(leaveI, gamma/(gamma+mu));
            enterS = leaveS - enterI + leaveI - enterR + leaveR;
            
            // Update state for next time step according to model
            SS(istep+1, jpatch) = iS1 + enterS - leaveS;  
            II(istep+1, jpatch) = iI1 + enterI - leaveI;
            RR(istep+1, jpatch) = iR1 + enterR - leaveR;
        }
        
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
    
    return ret;
};

// [[Rcpp::export]]
double seasonal_cosine(double t, List params) {
    double R0 = params["R0"];
    double gamma = params["gamma"];
    double pop = params["pop"];
    double b1 = params["b1"];
    double pi = 3.141592653589793;
    
    double b0 = R0 * gamma/pop;
    
    double ret = b0 * (1 + b1 * cos (2 * pi * t));
    
    return ret;
}

// [[Rcpp::export]]
double term_time(double t, List params) {
  double R0 = params["R0"];
  double gamma = params["gamma"];
  double pop = params["pop"];
  double b1 = params["b1"];
  // proportion of days in school year (weekends count, holidays do not)
  double p = 0.7589;
  
  // mean transmission rate
  double meantrans = R0 * gamma / pop;
  
  // seasonal amplitude
  double a = b1;
  
  double ret;
  
  // Reduce to a yearly loop
  // Data will repeat itself in simplified calendar
  
  t = (t - floor(t)) * 365;
  
  if(t >= 7 && t <= 68){
    ret = (1 + 2*(1-p)*a)*meantrans;
  } 
  else if(t >= 76 && t <= 182) {
    ret = (1 + 2*(1-p)*a)*meantrans;
  } 
  else if(t >= 246 && t <= 356) {
    ret = (1 + 2*(1-p)*a)*meantrans;
  } 
  else {
    ret = (1 - 2*p*a)*meantrans;
  }
  
  return ret;
  
}
