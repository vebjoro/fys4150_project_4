#include "state.hpp"
#include <armadillo>
#include <random>
#include <chrono>
#include <cmath>
#include <iostream>
#include <map>

// Constructor
State::State(int size, double temp, int seed)
{
  L = size;
  T = temp;
  uniform_dist = std::uniform_int_distribution<int>(1, L);
  uniform_real = std::uniform_real_distribution<>(0, 1);
  generator.seed(seed);
  M = 0; // total_magnetization() will set this
  E = 0; // total_energy() will set this
  E2 = 0;
  M2 = 0;
  dE = 0;
  dM = 0;

  // State matrix
  N = L * L;
  S = arma::mat(L + 2, L + 2);

  // Calculate five possible relative probabilities for acceptance
  for (int i = -2; i <= 2; i++)
  {
    prob_dict[i * 4.] = std::exp(-i * 4. / T);
  }
}

// Initialize random state
void State::init_random_state()
{

  S = arma::sign(S.randu() - 0.5);
  S.replace(0, 1);

  // Initialize periodic boundary conditions
  make_periodic();
}

// Initialize ordered positive state
void State::init_ordered_state()
{
  S = arma::ones(L + 2, L + 2);

  // Initialize periodic boundary conditions
  make_periodic();
}

void State::MC_cycle_sampling(int &j)
{
  // Assumes energy and magnetization of the previous microstate to be known.
  // Index j is the cycle number, the samples from each cycle is written to the
  // j-th element of vectors E_vec and M_vec.

  for (int i = 0; i < N; i++)
  {
    flip_random_spinn();
    E += dE;
    M += dM;
  }

  E_vec(j) = E;
  M_vec(j) = std::abs(M);
}

void State::MC_burn_in(int k)
{
  // Do N potential flips per MC-cycle, k is thecnumber of cycles.
  for (int i = 0; i < N * k; i++)
  {
    flip_random_spinn();
  }
}

void State::initialize_containers(int &n)

{
  // Initialize containers for energy and magnetization
  E_vec = arma::vec(n);
  M_vec = arma::vec(n);
}

void State::make_periodic()
{
  // Copy edges to ensure periodic boundary conditions
  arma::vec p_right = S.col(1);
  arma::vec p_left = S.col(L);
  S.col(0) = p_left;
  S.col(L + 1) = p_right;
  arma::rowvec p_top = S.row(1);
  arma::rowvec p_bottom = S.row(L);
  S.row(0) = p_bottom;
  S.row(L + 1) = p_top;

  // Set the corners to zero for visual purposes
  // S(0, 0) = 0;
  // S(L + 1, L + 1) = 0;
  // S(L + 1, 0) = 0;
  // S(0, L + 1) = 0;
}

// Flip a random spin according to the Boltzmann distribution for the given temperature
void State::flip_random_spinn()
{
  // Choose random spin
  int index_1 = uniform_dist(generator);
  int index_2 = uniform_dist(generator);

  // Difference in energy if spin is flipped
  double energy_diff = delta_E(index_1, index_2);

  // Difference in magnetization if spin is flipped
  double mag_diff = -2 * S(index_1, index_2);

  // Acceptance probability
  double rel_prob = prob_dict[energy_diff];

  // Accept or reject flip according to Boltzmann distribution
  double A = std::min(1., rel_prob);
  double r = uniform_real(generator);
  if (A > r)
  {
    S(index_1, index_2) *= -1;
    make_periodic();
    dE = energy_diff; // update change in system
    dM = mag_diff;
  }
  else
  {
    dE = 0;
    dM = 0;
  }
}

// Calculate the energy difference if a spin is flipped
double State::delta_E(int &index_1, int &index_2)
{
  double e_before = 0;
  e_before -= S(index_1, index_2) * S(index_1 + 1, index_2);
  e_before -= S(index_1, index_2) * S(index_1 - 1, index_2);
  e_before -= S(index_1, index_2) * S(index_1, index_2 + 1);
  e_before -= S(index_1, index_2) * S(index_1, index_2 - 1);

  return -2 * e_before;
}

// Calculate the total energy of one state
void State::total_energy()
{
  double energy = 0;
  for (int i = 1; i <= L; i++)
  {
    for (int j = 1; j <= L; j++)
    {
      energy -= S(i, j) * S(i + 1, j) + S(i, j) * S(i, j + 1);
    }
  }
  E = energy;
}

// Calculate the total magnetization of one state
void State::total_magnetization()
{
  M = arma::accu(S.rows(1, L).cols(1, L));
}

// Calculate the mean specific heat capacity
double State::specific_heat_capacity()
{
  double c = (1. / State::L / State::L) * (arma::mean(E_vec % E_vec) - arma::mean(E_vec) * arma::mean(E_vec));
  return c;
}

// Calculate the mean magnetic susceptibility
double State::magnetic_susceptibility()
{
  double x = (1. / State::L / State::L) * (arma::mean(M_vec % M_vec) - arma::mean(M_vec) * arma::mean(M_vec));
  return x;
}
