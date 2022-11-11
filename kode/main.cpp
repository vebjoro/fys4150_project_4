#include <iostream>
#include <armadillo>
#include "state.hpp"

int main(int argc, char *argv[])
{

  // Initialize 2x2 state
  int L = 2;
  double T = 1;
  int N = L * L;
  State state = State(L, T, 1);

  // Analytic solution 2x2
  double Z_1 = (2 * std::exp(8) + std::exp(-8) + 12);
  double analytic_energy = 4 * (std::exp(-8) - std::exp(8)) / Z_1;
  double analytic_magnetization = (8 * std::exp(8) + 16) / Z_1 / N;

  // Print analytic solution
  std::cout << "E_analytic: " << analytic_energy << std::endl;
  std::cout << "M_analytic: " << analytic_magnetization << std::endl;

  // Monte Carlo solution 2x2
  int n_cycles = 1000000;
  state.initialize_containers(n_cycles);

  for (int j = 0; j < n_cycles; j++)
  {
    state.MC_cycle_sampling(j);
  }

  // Print Monte Carlo solution
  double e_mean = arma::mean(state.e);
  double m_mean = arma::mean(state.m);
  std::cout << "--------------------------------------" << std::endl;
  std::cout << "E_calculated: " << e_mean << std::endl;
  std::cout << "M_calculated: " << m_mean << std::endl;

  return 0;
}
