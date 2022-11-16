#include <iostream>
#include <iomanip>
#include <armadillo>
#include "state.hpp"
#include "progressbar.hpp"

int main(int argc, char *argv[])
{
  std::string outfile;

  // Initialize 2x2 state
  int L = 2;
  double T = 1;
  int N = L * L;
  State state = State(L, T, 1);
  state.init_random_state();

  // Analytic solution 2x2
  double Z_1 = (2 * std::exp(8) + std::exp(-8) + 12);
  double e_analytic = 4 * (std::exp(-8) - std::exp(8)) / Z_1;
  double m_analytic = (8 * std::exp(8) + 16) / Z_1 / N;

  double expected_E = (-16 * std::exp(8) + 16 * std::exp(-8)) / Z_1;
  double expected_E2 = (128 * std::exp(8) + 128 * std::exp(-8)) / Z_1;

  double expected_M = (8 * std::exp(8) + 16) / Z_1;
  double expected_M2 = (32 * std::exp(8) + 32) / Z_1;

  double Cv_analytic = (1. / 4. * (expected_E2 - expected_E * expected_E)); // Specific heat capacity
  double X_analytic = (1. / 4. * (expected_M2 - expected_M * expected_M));  // Magnetic susceptibility

  // Monte Carlo solution 2x2
  int n_cycles = 100000;
  state.initialize_containers(n_cycles);
  state.MC_burn_in(1000);
  state.total_energy(); // calculating energy and mag of first microstate
  state.total_magnetization();

  for (int j = 0; j < n_cycles; j++)
  {
    state.MC_cycle_sampling(j);
  }

  // Print Monte Carlo solution
  double e_numeric = arma::mean(state.E_vec / N);
  double m_numeric = arma::mean(state.M_vec / N);
  double Cv_numeric = state.specific_heat_capacity();
  double X_numeric = state.magnetic_susceptibility();

  std::cout << "2x2 comparison" << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "      Analytic solution   Monte Carlo solution" << std::endl;
  std::cout << "e  : " << e_analytic << std::setw(21) << e_numeric << std::endl;
  std::cout << "m  :  " << m_analytic << std::setw(21) << m_numeric << std::endl;
  std::cout << "Cv :  " << Cv_analytic << std::setw(21) << Cv_numeric << std::endl;
  std::cout << "X  :  " << X_analytic << std::setw(21) << X_numeric << std::endl;

  return 0;
}
