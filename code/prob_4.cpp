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

  // Random initial state
  state.init_random_state();

  // Analytic solution 2x2
  double Z_1 = (2 * std::exp(8) + 2 * std::exp(-8) + 12);
  double e_analytic = 4 * (std::exp(-8) - std::exp(8)) / Z_1; // Expected value of energy per spin
  double m_analytic = (8 * std::exp(8) + 16) / Z_1 / N;       // Expected value of magnetization per spin

  double expected_E = (-16 * std::exp(8) + 16 * std::exp(-8)) / Z_1;
  double expected_E2 = (128 * std::exp(8) + 128 * std::exp(-8)) / Z_1;

  double expected_M = (8 * std::exp(8) + 16) / Z_1;
  double expected_M2 = (32 * std::exp(8) + 32) / Z_1;

  double Cv_analytic = (1. / 4. * (expected_E2 - expected_E * expected_E)); // Specific heat capacity
  double X_analytic = (1. / 4. * (expected_M2 - expected_M * expected_M));  // Magnetic susceptibility

  // Monte Carlo solution 2x2
  int n_cycles = 1000000;
  state.initialize_containers(n_cycles);
  state.MC_burn_in(10);
  state.total_energy(); // calculating energy and mag of first microstate
  state.total_magnetization();

  // Run Monte Carlo cycles
  for (int j = 0; j < n_cycles; j++)
  {
    state.MC_cycle_sampling(j);
  }

  // Calculate relative error
  double e_numeric = arma::mean(state.E_vec / N);
  double m_numeric = arma::mean(state.M_vec / N);
  double Cv_numeric = state.specific_heat_capacity();
  double X_numeric = state.magnetic_susceptibility();

  double e_relerr = std::abs((e_numeric - e_analytic)/e_analytic);
  double m_relerr = std::abs((m_numeric - m_analytic)/m_analytic);
  double Cv_relerr = std::abs((Cv_numeric - Cv_analytic)/Cv_analytic);
  double X_relerr = std::abs((X_numeric - X_analytic)/X_analytic);


  // Print Monte Carlo solution
  std::cout << "2x2 comparison" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "      Analytic solution   Monte Carlo solution   Relative error" << std::endl;
  std::cout << "e  : " << e_analytic << std::setw(21) << e_numeric << std::setw(21) << e_relerr << std::endl;
  std::cout << "m  :  " << m_analytic << std::setw(21) << m_numeric  << std::setw(20) << m_relerr << std::endl;
  std::cout << "Cv :  " << Cv_analytic << std::setw(21) << Cv_numeric  << std::setw(17) << Cv_relerr << std::endl;
  std::cout << "X  :  " << X_analytic << std::setw(21) << X_numeric  << std::setw(14) << X_relerr << std::endl;

  return 0;
}
