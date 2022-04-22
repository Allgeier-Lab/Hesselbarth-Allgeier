// [[Rcpp::depends(RcppProgress)]]

#include <Rcpp.h>
#include <chrono>
#include <random>
#include <progress.hpp>
#include <progress_bar.hpp>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector random_move(int n, int move_residence) {
  
  // progress bar
  Progress prgss(n, true);
  
  // init vector for results
  Rcpp::NumericVector result (n, 0.0);
  
  // loop through all repetitions
  for (int i = 0; i < n; i++) {
    
    // check if stopped
    if (Progress::check_abort()) Rcpp::stop("Stopped.");
    
    // loop through all time steps from 0 to maximum time
    for (int j = 0; j <= move_residence; j++) {
    
      // get random time-based seed 
      unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
      
      // set up random number generator
      std::mt19937 generator(seed);
      
      // random number between 0.0 and 1.0
      std::uniform_real_distribution<double> distribution(0.0, 1.0);
      
      // draw random number
      double prob_random = distribution(generator);
      
      // calculate movement probability
      double prob_move = j / (double) move_residence;
      
      // check if random number is bigger than move probs; nothing happens
      if (prob_random > prob_move) {
        
        continue;
      
      // this is the case individual would move
      } else {
        
        // save current time step for repetition
        result[i] = j;
        
        break;
      
      }
    }
    
    // increase progress bar
    prgss.increment();
    
  }

  // return result
  return result;
  
}
