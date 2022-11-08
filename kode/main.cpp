#include <iostream>
#include <armadillo>
#include "state.hpp"

int main(int argc, char *argv[]) {

  int L = 2;
  double T = 1;
  int N = L*L;
  State state = State(L, T, 1);
  std::cout << state.S << std::endl;


  int n_cycles = 10;
  state.initialize_containers(n_cycles);
  for (int j; j < n_cycles; j++){
    state.MC_cycle_sampling(j);
    std::cout << state.chai[j] << std::endl;
  }

  //std::cout << state.S << std::endl;
  return 0;
}
