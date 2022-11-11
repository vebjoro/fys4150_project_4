#ifndef STATE_H_
#define STATE_H_

#include <armadillo>

// Structure
struct State
{

  int L;
  int N;
  double T;
  arma::mat S;
  int seed;
  std::mt19937 generator;
  std::uniform_int_distribution<int> uniform_dist;
  std::uniform_real_distribution<> uniform_real;

  double E;
  double M;
  double E2;
  double M2;

  arma::vec e;
  arma::vec m;
  arma::vec Cv;
  arma::vec chi;

  State(int L, double temp, int seed);

  void init_random_state();

  void init_uniform_state();

  void MC_cycle_sampling(int &j);

  void initialize_containers(int &n);

  double flip_random_spinn();

  void make_periodic();

  double delta_E(int &index_1, int &index_2);

  double total_energy();

  double total_magnetization();

  double specific_heat_capacity();

  double magnetic_susceptibility();
};

#endif // STATE_H_
