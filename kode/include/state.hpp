#ifndef STATE_H_
#define STATE_H_

#include <armadillo>

// Structure
struct State{

  int L ;
 double T;
  arma::mat S;
  int seed;
  std::mt19937 generator;
  std::uniform_int_distribution<int> uniform_dist;
  std::uniform_real_distribution<> uniform_real;

  State(int L, double temp, int seed);

  void flip_random_spinn();

  void make_periodic();

  double delta_E(int &index_1, int &index_2);


};


#endif // STATE_H_
