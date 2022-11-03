#ifndef STATE_H_
#define STATE_H_

#include <armadillo>

// Structure
struct State{

  int L ;
 double T;
  arma::mat S;


  State(int L, double temp);

  void flip_random_spinn();

  double delta_E(int &index_1, int &index_2);



};


#endif // STATE_H_
