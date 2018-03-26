#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List SIRstochastic(List params) {
  
  // Number of time steps
  int nsteps = as<int>(params["nsteps"]);
  
  // Number of patches
  int mpatches = as<int>(params["mpatches"]);
  
  // Pull list of initial states
  List init = params["init"];
  
  // Initalize matrices to hold results (time steps as rows, columns as patches) - makes matrix of zeroes
  NumericMatrix SS(nsteps, mpatches);
  NumericMatrix II(nsteps, mpatches);
  NumericMatrix RR(nsteps, mpatches);
  
  // Get the initial state values
  double Sinit = init["S"];
  double Iinit = init["I"];
  double Rinit = init["R"];
  
  // Set initial conditions in the first row
  for (int i = 0; i < mpatches; i++) {
    SS(0, i) = Sinit;
    II(0, i) = Iinit;
    RR(0, i) = Rinit;
  }
  
  // Initialize time vector
  NumericVector time(nsteps);
  
  // Pull out parameter values
  double dt = params["dt"]; 
  double mu = params["mu"];
  double beta = params["beta"];
  double gamma = params["gamma"];
  double pop = params["pop"];
  double births = mu*pop*dt;
  double recoveryrate = 1 - exp(-(gamma+mu)*dt);
  double deathrate = 1 - exp(-mu*dt);
  
  // Connectivity matrix
  NumericMatrix m = params["conmat"];
  
  // Fill connectivity matrix (currently with 0.5)
  for (int i = 0; i < mpatches*mpatches; i++) {
    m[i] = 0.5;
  }
  
  // Iterate over time steps
  for (int istep = 0; istep < (nsteps-1); istep++) {
    
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
      
      // The proportion that leaves the susceptible class
      double S_leave = R::rbinom(iS1, 1 - exp(-(beta*transsum + mu)*dt));
      
      // The proportion that leaves the infected class
      double I_leave = R::rbinom(iI1, recoveryrate);
      
      // The proportion that leaves the recovered class
      double R_leave = R::rbinom(iR1, deathrate);
      
      // Number of new infected individuals produced in time interval (t, t+dt)
      double i_k = R::rbinom(S_leave, ((beta*transsum)/(beta*transsum + mu)));
      
      // Number of new recovered individuals produced in time interval (t, t+dt)
      double r_k = R::rbinom(R_leave, (gamma/(gamma+mu)));
      
      // Number of new susceptible individuals produced in time interval (t, t+dt)
      double b_k = S_leave - i_k + I_leave - r_k;
      
      // Update state for next time step according to model
      SS(istep+1, jpatch) = iS1 + b_k - S_leave;  
      II(istep+1, jpatch) = iI1 + i_k - I_leave;
      RR(istep+1, jpatch) = iR1 + r_k - R_leave;
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