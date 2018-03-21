#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List SIRmodel(List params) {
  
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
    
    // The proportion that leave the susceptible class
    double transmission = 1- exp((-beta*transsum -mu)*dt);
    
    //i_i(t) - the number of individuals going from S to I
    double transratio = beta*transsum*transmission*iS1/(beta*transsum + mu);
    
    // r_i(t) - the number of individuals going from I to R
    double recoveryratio = gamma*recoveryrate*iI1/(gamma+mu); 
      
    // Update state for next time step according to model
    SS(istep+1, jpatch) = iS1 + births - transmission*iS1;  
    II(istep+1, jpatch) = iI1 + transratio - recoveryrate*iI1;
    RR(istep+1, jpatch) = iR1 + recoveryratio - deathrate*iR1;
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