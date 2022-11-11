#include <iostream>
#include <armadillo>
#include "state.hpp"

int main(int argc, char *argv[]) {

  int L = 2;
  double T = 1;
  int N = L*L;
  State state = State(L, T, 1);
  //std::cout << state.S << std::endl;

  double Z_1 = (2*std::exp(8) + std::exp(-8) + 12);
  double analytic_energy = 4*(std::exp(-8) - std::exp(8))/Z_1;
  double analytic_magnetization = (8*std::exp(8) + 16)/Z_1/N;
  std::cout << "E_analytic: " << analytic_energy << std::endl;
  std::cout << "M_analytic: " << analytic_magnetization << std::endl;

  int n_cycles = 1000000;
 state.initialize_containers(n_cycles);
  for (int j; j < n_cycles; j++){
 //   std::cout << "New cycle" << std::endl;
    state.MC_cycle_sampling(j);
  }

  double e_mean = arma::mean(state.e);
  double m_mean = arma::mean(state.m);
  std::cout << "--------------------------------------" << std::endl;
  //std::cout << state.e << std::endl;
 std::cout << "E_calculated: " << e_mean << std::endl;
  //std::cout << state.m << std::endl;
  std::cout << "M_calculated: " <<  m_mean << std::endl;

  //std::cout << state.S << std::endl;
  return 0;
}
