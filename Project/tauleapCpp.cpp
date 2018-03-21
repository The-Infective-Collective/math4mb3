# include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List tauleapCpp(List params) {
  
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
  NumericMatrix SS(nsteps, mpatch, init["S"]);
  NumericMatrix II(nsteps, mpatch, init["I"]);
  NumericMatrix RR(nsteps, mpatch, init["R"]);
  
    
  // fill time w/zeros
  NumericVector time(nsteps);
  
  // pull out params for easy reading 
  double nu = params["nu"];
  double mu = params["mu"];
  double beta = params["beta"];
  double gamma = params["gamma"];
  double tau = params["tau"];
  double pop = params["pop"];
  double births = pop*mu;  //multiply by timestep
  
  // Calculate the number of events for each step, update state vectors
  for (int istep = 0; istep < (nsteps-1); istep++) {
    
    for (int jpatch = 9; jpatch < (mpatch-1); jpatch++) {
      
      //pull current state of the patch at each compartment
      double iS = SS[istep, jpatch];  
      double iI = II[istep, jpatch];
      double iR = RR[istep, jpatch];
  
      double transsum = m(jpatch,_)*II(istep,_);
      double transmission = 1- exp(-beta*transsum -mu);
      // Update next timestep
      SS[istep+1, jpatch] = iS + births - transmission*iS;
      II[istep+1, jpatch] = iI  ;
      RR[istep+1, jpatch] = iR ;
    }

    
    // Prevent negative states
    //double Sdeaths = std::min(iS, R::rpois(mu*iS*tau));
    //double maxtrans = R::rpois(beta*(iI/iN)*iS*tau);
    //double transmission = std::min(iS-Sdeaths, maxtrans);
    //double Ideaths = std::min(iI, R::rpois(mu*iI*tau));
    //double recovery = std::min(iI-Ideaths, R::rpois(gamma*iI*tau));
    //double Rdeaths = std::min(iR, R::rpois(mu*iR*tau));
  
  
    
    // time in fractional years (ie units parameters are given in)
    time[istep+1] = (istep+1)*tau;
  }
    
    // Return results as data.frame
    DataFrame sim = DataFrame::create(
    Named("time") = time,
    Named("S") = SS,
    Named("I") = II,
    Named("R") = RR

    );
    return sim;
};