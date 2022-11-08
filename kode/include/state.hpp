#ifndef STATE_H_
#define STATE_H_

#include <armadillo>

// Structure
struct State{

  int L ;
 int N;
  double T;
  arma::mat S;
  int seed;
  std::mt19937 generator;
  std::uniform_int_distribution<int> uniform_dist;
  std::uniform_real_distribution<> uniform_real;

  double E_sample;
  double M_sample;
 double E2_sample;
  double M2_sample;

  arma::vec e;
  arma::vec m;
  arma::vec Cv;
  arma::vec chai;

  State(int L, double temp, int seed);

  void MC_cycle_sampling(int &j);

  void initialize_containers(int &n);

  void flip_random_spinn();

  void make_periodic();

  double delta_E(int &index_1, int &index_2);

  double total_energy();

  double total_magnetization();

  double specific_heat_capacity();

  double magnetic_susceptibility();



};


#endif // STATE_H_
