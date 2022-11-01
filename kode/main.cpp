#include <iostream>
#include <armadillo>
#include "state.hpp"

int main(int argc, char *argv[]) {


  State state = State(2, 1);
  std::cout << state.S << std::endl;


  // int n_cycles = 10;
  // for (int j; j < n_cycles; j++){
  //   for (int i = 0; i<N; i++){
  //       flip_random_spinn(state, L, T);
  //       //std::cout << state << std::endl;
  //    }
  // }

  return 0;
}
