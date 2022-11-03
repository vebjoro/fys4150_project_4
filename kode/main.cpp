#include <iostream>
#include <armadillo>
#include "state.hpp"

int main(int argc, char *argv[]) {

  int L = 2;
  double T = 1;
  int N = 100;
  State state = State(L, T);
  std::cout << state.S << std::endl;


  // int n_cycles = 10;
  //for (int j; j < n_cycles; j++){
    for (int i = 0; i<N; i++){
        state.flip_random_spinn();
        //std::cout << state << std::endl;
     }
  // }

  std::cout << state.S << std::endl;
  return 0;
}
