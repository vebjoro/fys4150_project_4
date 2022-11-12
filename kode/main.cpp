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
  state.init_random_state();

  // Analytic solution 2x2
  double Z_1 = (2 * std::exp(8) + std::exp(-8) + 12);
  double analytic_energy = 4 * (std::exp(-8) - std::exp(8)) / Z_1;
  double analytic_magnetization = (8 * std::exp(8) + 16) / Z_1 / N;

  // Print analytic solution
 std::cout << "E_analytic: " << analytic_energy << std::endl;
 std::cout << "M_analytic: " << analytic_magnetization << std::endl;

  // Monte Carlo solution 2x2
  int n_cycles = 100000;
  state.initialize_containers(n_cycles);
  state.MC_burn_in(1000);
  state.total_energy(); //calculating energy and mag of first microstate
  state.total_magnetization();

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

//   // Initialize 20x20 state for T = 1 J/k_B and T = 2.4 J/k_B
  L = 20;
  n_cycles = 10000;

  T = 1;
  State state_20_t01_Random = State(L, T, 1);
  state_20_t01_Random.init_random_state(); // Random initial state
  state_20_t01_Random.initialize_containers(n_cycles);
   state_20_t01_Random.total_energy();
   state_20_t01_Random.total_magnetization();


  State state_20_t01_Ordered = State(L, T, 1);
  state_20_t01_Ordered.init_ordered_state(); // Ordered initial state
  state_20_t01_Ordered.initialize_containers(n_cycles);
  state_20_t01_Ordered.total_energy();
   state_20_t01_Ordered.total_magnetization();

  T = 2.4;
  State state_20_t24_Random = State(L, T, 1);
  state_20_t24_Random.init_random_state(); // Random initial state
  state_20_t24_Random.initialize_containers(n_cycles);
  state_20_t24_Random.total_energy();
   state_20_t24_Random.total_magnetization();


  State state_20_t24_Ordered = State(L, T, 1);
  state_20_t24_Ordered.init_ordered_state(); // Ordered initial state
  state_20_t24_Ordered.initialize_containers(n_cycles);
   state_20_t24_Ordered.total_energy();
   state_20_t24_Ordered.total_magnetization();


  // Simulate states 20x20
  arma::mat out_mean_E = arma::zeros(4, n_cycles); // 4 states, n_cycles (Energy)
  arma::mat out_mean_M = arma::zeros(4, n_cycles); // 4 states, n_cycles (Magnetization)
  arma::mat out_pdf_E = arma::zeros(4, n_cycles); // 4 states, n_cycles (Energy)
  arma::mat out_pdf_M = arma::zeros(4, n_cycles);


  //sample expectation value to check convergence
  for (int j = 0; j < n_cycles; j++)
  {
    state_20_t01_Random.MC_cycle_sampling(j);
    state_20_t01_Ordered.MC_cycle_sampling(j);
    state_20_t24_Random.MC_cycle_sampling(j);
     state_20_t24_Ordered.MC_cycle_sampling(j);

    out_mean_E(0, j) = arma::mean(state_20_t01_Random.e(arma::span(0, j)));
    out_mean_E(1, j) = arma::mean(state_20_t01_Ordered.e(arma::span(0, j)));
    out_mean_E(2, j) = arma::mean(state_20_t24_Random.e(arma::span(0, j)));
    out_mean_E(3, j) = arma::mean(state_20_t24_Ordered.e(arma::span(0, j)));

    out_mean_M(0, j) = arma::mean(state_20_t01_Random.m(arma::span(0, j)));
    out_mean_M(1, j) = arma::mean(state_20_t01_Ordered.m(arma::span(0, j)));
    out_mean_M(2, j) = arma::mean(state_20_t24_Random.m(arma::span(0, j)));
    out_mean_M(3, j) = arma::mean(state_20_t24_Ordered.m(arma::span(0, j)));
  }


  // Save to file
  std::string outfile;

  outfile = "plot/binary_data/20x20_E_mean.bin";
  out_mean_E.save(outfile, arma::arma_binary);

  outfile = "plot/binary_data/20x20_M_mean.bin";
  out_mean_M.save(outfile, arma::arma_binary);

//sample pdf, assume 20x20 states to be burned in

   return 0;
 }
