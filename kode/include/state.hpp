#ifndef STATE_H_
#define STATE_H_

#include <armadillo>
#include <map>

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

  std::map<double, double> prob_dict;

  double E;
  double M;
  double E2;
  double M2;
  double dE;
  double dM;

  arma::vec e;
  arma::vec m;
  arma::vec E_vec;
  arma::vec M_vec;
  arma::vec Cv;
  arma::vec chi;

  State(int L, double temp, int seed);

  void init_random_state();

  void init_ordered_state();

  void MC_cycle_sampling(int &j);

  void MC_burn_in(int k);

  void initialize_containers(int &n);

  void flip_random_spinn();

  void make_periodic();

  double delta_E(int &index_1, int &index_2);

  void total_energy();

  void total_magnetization();

  double specific_heat_capacity();

  double magnetic_susceptibility();
};

#endif // STATE_H_
