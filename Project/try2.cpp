#include <Rcpp.h>
using namespace Rcpp;

// This function will be used in R! Evaluates the number of events 
// and updates the states at each time step
//
// [[Rcpp::export]]
List SIRmodel(List params) {
  
  // chained operations are tricky in cpp
  // pull out list w/in list into its own object
  List init = params["init"];
  
  // use Rcpp as() function to "cast" R vector to cpp scalar
  int nsteps = as<int>(params["nsteps"]);
  int mpatch = as<int>(params["mpatch"]);
  
  // initialize each state vector in its own vector
  // set all vals to initial vals
  //
  // I use doubles (NumericVector) rather than 
  // ints (IntegerVector), since rpois returns double,
  // and the domain of double is a superset of int
  NumericVector SS(nsteps, init["S"]);
  NumericVector II(nsteps, init["I"]);
  NumericVector RR(nsteps, init["R"]);
  
  
  NumericMatrix SS1(nsteps, mpatch);
  NumericMatrix II1(nsteps, mpatch);
  NumericMatrix RR1(nsteps, mpatch);
  
  double Sinit = init["S"];
  double Iinit = init["I"];
  double Rinit = init["R"];
  int size = SS1.nrow()*SS1.ncol();
  
  for (int i = 0; i < size; i++) {
    SS1[i] = Sinit;
    II1[i] = Iinit;
    RR1[i] = Rinit;
  }
  
  // fill time w/zeros
  NumericVector time(nsteps);
  
  // pull out params for easy reading 
  double dt = params["dt"];
  double mu = params["mu"];
  double beta = params["beta"];
  double gamma = params["gamma"];
  double pop = params["pop"];
  
  double births = mu*pop*dt;
  double recoveryrate = 1 - exp(-(gamma+mu)*dt);
  double deathrate = 1 - exp(-mu*dt);
  
  NumericMatrix m(mpatch, mpatch);
  
  for (int i = 0; i < mpatch*mpatch; i++) {
    m[i] = 0.5;
  }
  
  // Calculate the number of events for each step, update state vectors
  for (int istep = 0; istep < (nsteps-1); istep++) {
    
    for (int jpatch = 0; jpatch < mpatch; jpatch++) {
      
      double iS1 = SS1(istep, jpatch);  
      double iI1 = II1(istep, jpatch);
      double iR1 = RR1(istep, jpatch);
      
      // perform dot product for transmission rate
      double transsum = 0;
      
      for (int i=0; i < (mpatch); i++) {
        transsum += m(jpatch, i)*II1(istep,i);
      }
    
    double transmission = 1- exp((-beta*transsum -mu)*dt);
    
    //i_i(t)
    double transratio = beta*transsum*transmission*iS1/(beta*transsum + mu);
    
    // r_i(t)
    double recoveryratio = gamma*recoveryrate*iI1/(gamma+mu); 
      
    SS1(istep+1, jpatch) = iS1 + births - transmission*iS1;  
    II1(istep+1, jpatch) = iI1 + transratio - recoveryrate*iI1;
    RR1(istep+1, jpatch) = iR1 + recoveryratio - deathrate*iR1;
    }
    
    // time in fractional years (ie units parameters are given in)
    time[istep+1] = (istep+1)*dt;
  }

  List ret;
  ret["time"] = time;
  ret["S"] = SS1;
  ret["I"] = II1;
  ret["R"] = RR1;
  ret["m"] = m;
  return ret;
};