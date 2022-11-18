#include <iostream>
#include <iomanip>
#include <armadillo>
#include "omp.h"
#include "state.hpp"
#include "progressbar.hpp"
#include <stdio.h>
#include <chrono>

int main(int argc, char *argv[])
{
  std::string outfile;

  // Phase transition investigation with OpenMP
  int number_of_temperatures = 12; // TODO: Discuss this
  int n_cycles = 1000000;          // TODO: Discuss this
  int n_burn = 2500;
  // arma::vec T_vec = arma::linspace(2.1, 2.4, number_of_temperatures); // Temperatures to be investigated
  arma::vec T_vec = arma::linspace(2.26, 2.34, number_of_temperatures); // Temperatures to be investigated

  arma::mat e_out = arma::zeros(4, number_of_temperatures); // 4 lattice sizes
  arma::mat m_out = arma::zeros(4, number_of_temperatures);
  arma::mat Cv_out = arma::zeros(4, number_of_temperatures);
  arma::mat X_out = arma::zeros(4, number_of_temperatures);

  // timing
  auto begin = std::chrono::high_resolution_clock::now();

#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < number_of_temperatures; i++)
    {

      State state_40 = State(40, T_vec(i), 1);
      state_40.init_random_state(); // Random initial state
      state_40.MC_burn_in(n_burn);  // TODO: Discuss this
      state_40.initialize_containers(n_cycles);
      state_40.total_energy();
      state_40.total_magnetization();

      State state_60 = State(60, T_vec(i), 1);
      state_60.init_random_state(); // Random initial state
      state_60.MC_burn_in(n_burn);
      state_60.initialize_containers(n_cycles);
      state_60.total_energy();
      state_60.total_magnetization();

      State state_80 = State(80, T_vec(i), 1);
      state_80.init_random_state(); // Random initial state
      state_80.MC_burn_in(n_burn);
      state_80.initialize_containers(n_cycles);
      state_80.total_energy();
      state_80.total_magnetization();

      State state_100 = State(100, T_vec(i), 1);
      state_100.init_random_state(); // Random initial state
      state_100.MC_burn_in(n_burn);
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

  auto end = std::chrono::high_resolution_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
  std::cout << "N_cycles: " << n_cycles << ", "
            << "T_steps: " << number_of_temperatures << std::endl;
  std::cout << "Elapsed time: " << elapsed.count() * 1e-9 << "s" << std::endl;

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
