#include <iostream>
#include <iomanip>
#include <armadillo>
#include "state.hpp"
#include "progressbar.hpp"

int main(int argc, char *argv[])
{
  std::string outfile;

  // Initialize 20x20 state for T = 1 J/k_B and T = 2.4 J/k_B
  int L = 20;
  int N = L * L;
  int n_cycles = 10000;

  double T = 1;
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

  // prob 5
  //  Simulate states 20x20
  arma::mat out_mean_E = arma::zeros(4, n_cycles); // 4 states, n_cycles (Energy)
  arma::mat out_mean_M = arma::zeros(4, n_cycles); // 4 states, n_cycles (Magnetization)
  arma::mat out_pdf_E = arma::zeros(4, n_cycles);  // 4 states, n_cycles (Energy)
  arma::mat out_pdf_M = arma::zeros(4, n_cycles);

  for (int j = 0; j < n_cycles; j++)
  {
    state_20_t01_Random.MC_cycle_sampling(j);
    state_20_t01_Ordered.MC_cycle_sampling(j);
    state_20_t24_Random.MC_cycle_sampling(j);
    state_20_t24_Ordered.MC_cycle_sampling(j);

    out_mean_E(0, j) = arma::mean(state_20_t01_Random.E_vec(arma::span(0, j)) / N);
    out_mean_E(1, j) = arma::mean(state_20_t01_Ordered.E_vec(arma::span(0, j)) / N);
    out_mean_E(2, j) = arma::mean(state_20_t24_Random.E_vec(arma::span(0, j)) / N);
    out_mean_E(3, j) = arma::mean(state_20_t24_Ordered.E_vec(arma::span(0, j)) / N);

    out_mean_M(0, j) = arma::mean(state_20_t01_Random.M_vec(arma::span(0, j)) / N);
    out_mean_M(1, j) = arma::mean(state_20_t01_Ordered.M_vec(arma::span(0, j)) / N);
    out_mean_M(2, j) = arma::mean(state_20_t24_Random.M_vec(arma::span(0, j)) / N);
    out_mean_M(3, j) = arma::mean(state_20_t24_Ordered.M_vec(arma::span(0, j)) / N);
  }

  // Save to file

  outfile = "plot/binary_data/20x20_E_mean.bin";
  out_mean_E.save(outfile, arma::arma_binary);

  outfile = "plot/binary_data/20x20_M_mean.bin";
  out_mean_M.save(outfile, arma::arma_binary);

  // prob 6
  //  sample pdf, assume 20x20 states to be burned in
  n_cycles = 100000;
  state_20_t01_Random.initialize_containers(n_cycles);
  state_20_t24_Random.initialize_containers(n_cycles);
  state_20_t01_Random.total_energy();
  state_20_t24_Random.total_energy();
  state_20_t01_Random.total_magnetization();
  state_20_t24_Random.total_magnetization();

  progressbar bar(n_cycles);
  std::cout << "Sampling pdf" << std::endl;
  for (int j = 0; j < n_cycles; j++)
  {
    state_20_t01_Random.MC_cycle_sampling(j);
    state_20_t24_Random.MC_cycle_sampling(j);
    bar.update();
  }

  outfile = "plot/binary_data/20x20_1_e.bin";
  arma::vec e = state_20_t01_Random.E_vec / N;
  e.save(outfile, arma::arma_binary);

  outfile = "plot/binary_data/20x20_24_e.bin";
  e = state_20_t24_Random.E_vec / N;
  e.save(outfile, arma::arma_binary);

  return 0;
}
