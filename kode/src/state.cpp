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
  // PDF
  uniform_dist = std::uniform_int_distribution<int>(1, L);
  uniform_real = std::uniform_real_distribution<>(0, 1);
  generator.seed(seed);
  M = 0; //total_magnetization();
  E = 0; //total_energy();
  E2 = 0;
  M2 = 0;
  dE = 0;
  dM = 0;

  // State matrix
  N = L * L;
  S = arma::mat(L + 2, L + 2);

  //calculate five possible relative probabilities for acceptance
  for (int i=-2; i<=2; i++){
    prob_dict[i*4.] = std::exp(-i*4. / T);
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

void State::MC_cycle_sampling(int &j){
  //assumes energy and magnetization of the previous microstate to be known
  //j is the cycle number, the samples from each cycle is written
  //to the j-th element of vectors e, m, Cv, chi, representing mean energy per spin, mean magnetization per
  //spin, specific heat capacity and magnetic susceptibility
 //total_energy();
 //total_magnetization();
 // std::cout << "-------------------------------------------" << std::endl;
  //std::cout << "Før: " << E << std::endl;
  for (int i = 0; i < N; i++)
  {
    flip_random_spinn();
    E +=dE;
    M+=dM;
  }
  //std::cout << "E+dE: " << E << std::endl;
  //total_energy();
  //std::cout <<"Total energi: "<< E << std::endl;
  //total_magnetization();

  E2 = std::pow(E, 2);
  M2 = std::pow(M, 2);

  // Mean over number of spins
  e(j) = E / N;
  m(j) = std::abs(M / N);

  Cv(j) = (specific_heat_capacity());
  chi(j) = (magnetic_susceptibility());
}


void State::MC_burn_in(int k){
  //do N potential flips per MC-cycle, k is the
  //number of cycles
  for (int i = 0; i < N*k; i++)
  {
    flip_random_spinn();
  }
}

void State::initialize_containers(int &n)
{
  e = arma::vec(n);
  m = arma::vec(n);
  Cv = arma::vec(n);
  chi = arma::vec(n);
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
  S(0, 0) = 0;
  S(L + 1, L + 1) = 0;
  S(L + 1, 0) = 0;
  S(0, L + 1) = 0;
}

// Flip a random spin according to the Boltzmann dist. for the given temperature
void State::flip_random_spinn()
{
  // Chose random spin
  int index_1 = uniform_dist(generator); // korleis kalle på arma-element med ein index?
  int index_2 = uniform_dist(generator);

  // std::cout << "Indeksar: " << index_1  << ", " << index_2 << std::endl;

  //energiforskjel om spinnen flippast
  double energy_diff = delta_E(index_1, index_2);
  double mag_diff = -2*S(index_1, index_2);
  //std::cout << mag_diff << std::endl;
  // std::cout << energy_diff << std::endl;
  double rel_prob =  prob_dict[energy_diff];
  //double rel_prob = std::exp(-energy_diff/T);
  //std::cout << "------------------------------------------------------------------------------" << std::endl;
   //std::cout << S.rows(index_1-1, index_1+1).cols(index_2-1, index_2+1) << std::endl;
  // Accept or reject
   double A = std::min(1., rel_prob);
   double r = uniform_real(generator);
    if (A > r){
      S(index_1, index_2) *= -1;
      make_periodic();
      dE = energy_diff; //update change in system
      dM = mag_diff;
    }
    else{
    dE = 0;
    dM = 0;
    }
    //std::cout << S.rows(index_1-1, index_1+1).cols(index_2-1, index_2+1) << std::endl;
    //std::cout << dE << std::endl;
}

// Calculate the energy difference if a spin is flipped
double State::delta_E(int &index_1, int &index_2)
{
  // std::cout << "Calculate energy difference" <<std::endl;
  double e_before = 0;
  e_before -= S(index_1, index_2) * S(index_1 + 1, index_2);
  e_before -= S(index_1, index_2) * S(index_1 - 1, index_2);
  e_before -= S(index_1, index_2) * S(index_1, index_2 + 1);
  e_before -= S(index_1, index_2) * S(index_1, index_2 - 1);

  return -2 * e_before;
}

// Calculate the total energy of one state
// store to member variable E
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
// store to member var M
void State::total_magnetization()
{
  M = arma::accu(S.rows(1, L).cols(1, L));
}

double State::specific_heat_capacity()
// Calculates the specific heat capacity for a state using samples from one MC cycle
{
  double c = 1 / (N * T * T) * (E2 / N - std::pow(E / N, 2));
  return c;
}

double State::magnetic_susceptibility()
// Calculates the magnetic susceptibility for a state using samples from one MC cycle
{
  double x = 1 / (N * T) * (M2 / N - std::pow(M / N, 2));
  return x;
}
