#include <iostream>
#include <armadillo>
#include "omp.h"
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

  // Initialize 20x20 state for T = 1 J/k_B and T = 2.4 J/k_B
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

  // sample pdf, assume 20x20 states to be burned in
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
  state_20_t01_Random.e.save(outfile, arma::arma_binary);

  outfile = "plot/binary_data/20x20_24_e.bin";
  state_20_t24_Random.e.save(outfile, arma::arma_binary);

  // Phase transition investigation with OpenMP
  int number_of_temperatures = 80;                                    // TODO: Discuss this
  n_cycles = 3000;                                                    // TODO: Discuss this
  arma::vec T_vec = arma::linspace(2.1, 2.4, number_of_temperatures); // Temperatures to be investigated

  arma::mat e_out = arma::zeros(4, number_of_temperatures); // 4 lattice sizes
  arma::mat m_out = arma::zeros(4, number_of_temperatures);
  arma::mat Cv_out = arma::zeros(4, number_of_temperatures);
  arma::mat X_out = arma::zeros(4, number_of_temperatures);

#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < number_of_temperatures; i++)

    {

      State state_40 = State(40, T_vec(i), 1);
      state_40.init_random_state(); // Random initial state
      state_40.initialize_containers(n_cycles);
      state_40.total_energy();
      state_40.total_magnetization();

      State state_60 = State(60, T_vec(i), 1);
      state_60.init_random_state(); // Random initial state
      state_60.initialize_containers(n_cycles);
      state_60.total_energy();
      state_60.total_magnetization();

      State state_80 = State(80, T_vec(i), 1);
      state_80.init_random_state(); // Random initial state
      state_80.initialize_containers(n_cycles);
      state_80.total_energy();
      state_80.total_magnetization();

      State state_100 = State(100, T_vec(i), 1);
      state_100.init_random_state(); // Random initial state
      state_100.initialize_containers(n_cycles);
      state_100.total_energy();
      state_100.total_magnetization();

      for (int j = 0; j < n_cycles; j++)
      {
        state_40.MC_cycle_sampling(j);
        state_60.MC_cycle_sampling(j);
        state_80.MC_cycle_sampling(j);
        state_100.MC_cycle_sampling(j);
      }

      e_out(0, i) = arma::mean(state_40.E_vec / 40 / 40);
      e_out(1, i) = arma::mean(state_60.E_vec / 60 / 60);
      e_out(2, i) = arma::mean(state_80.E_vec / 80 / 80);
      e_out(3, i) = arma::mean(state_100.E_vec / 100 / 100);

      m_out(0, i) = arma::mean(state_40.M_vec / 40 / 40);
      m_out(1, i) = arma::mean(state_60.M_vec / 60 / 60);
      m_out(2, i) = arma::mean(state_80.M_vec / 80 / 80);
      m_out(3, i) = arma::mean(state_100.M_vec / 100 / 100);

      Cv_out(0, i) = state_40.specific_heat_capacity();
      Cv_out(1, i) = state_60.specific_heat_capacity();
      Cv_out(2, i) = state_80.specific_heat_capacity();
      Cv_out(3, i) = state_100.specific_heat_capacity();

      X_out(0, i) = state_40.magnetic_susceptibility();
      X_out(1, i) = state_60.magnetic_susceptibility();
      X_out(2, i) = state_80.magnetic_susceptibility();
      X_out(3, i) = state_100.magnetic_susceptibility();
    }
  } // end pragma omp parallel

  outfile = "plot/binary_data/OMP_T_out.bin";
  T_vec.save(outfile, arma::arma_binary);

  outfile = "plot/binary_data/OMP_e_out.bin";
  e_out.save(outfile, arma::arma_binary);

  outfile = "plot/binary_data/OMP_m_out.bin";
  m_out.save(outfile, arma::arma_binary);

  outfile = "plot/binary_data/OMP_Cv_out.bin";
  Cv_out.save(outfile, arma::arma_binary);

  outfile = "plot/binary_data/OMP_X_out.bin";
  X_out.save(outfile, arma::arma_binary);

  return 0;
}