#include "state.hpp"
#include <armadillo>
#include <random>
#include <chrono>
#include <cmath>


//konstruktør
State::State(int size, double temp, int seed){
  L = size;
  T = temp;
  E_sample = 0;
  M_sample = 0;
  E2_sample = 0;
  M2_sample = 0;

  //pdf
 uniform_dist = std::uniform_int_distribution<int>(1, L);
  uniform_real = std::uniform_real_distribution<> (0, 1);
  generator.seed(seed);

  //tilstandsmatrise
   N = L*L;
  S = arma::mat(L+2, L+2);
   S = arma::sign(S.randu()-0.3);
  S.replace(0, 1);

  //periodisitet
  make_periodic();
}



void State::MC_cycle_sampling(int &j)
{
  E_sample = 0;
  M_sample = 0;
  for (int i = 0; i<N; i++){
        flip_random_spinn();
        E_sample += total_energy();
        M_sample += total_magnetization();
    }
  e(j) = (E_sample/std::pow(N, 2)); //mean over number of spins and number of (potential) flips
  m(j) = (M_sample/std::pow(N, 2));
  Cv(j) = (specific_heat_capacity());
  chai(j) = (magnetic_susceptibility());
}

void State::initialize_containers(int &n)
{
  e = arma::vec(n);
  m = arma::vec(n);
  Cv = arma::vec(n);
  chai = arma::vec(n);
}


void State::make_periodic()
  {
    //kopierar rendene til tilstanden for å forsikre periodisitet

  arma::vec p_right = S.col(1);
  arma::vec p_left = S.col(L);
  S.col(0) = p_left;
  S.col(L+1) = p_right;
 arma::rowvec p_top = S.row(1);
  arma::rowvec p_bottom = S.row(L);
  S.row(0) = p_bottom;
  S.row(L+1) = p_top;
  //lettare å visualisere
  S(0, 0) = 0;
   S(L+1, L+1) = 0;
   S(L+1, 0) = 0;
   S(0, L+1) = 0;
}


//flippar ein spinn
void State::flip_random_spinn()
{
  //std::cout << "Calling flip spinn" << std::endl;

  //brukar generatoren til å velge ein tilfeldig spinn
  int index_1  = uniform_dist(generator); // korleis kalle på arma-element med ein index?
  int index_2  = uniform_dist(generator);
  //std::cout << "Indeksar: " << index_1  << ", " << index_2 << std::endl;

  //energiforskjel om spinnen flippast
   double energy_diff = delta_E(index_1, index_2);
   std::cout << energy_diff << std::endl;
   double rel_prob = std::exp(-energy_diff/T); //double check this!

   //godkjenningssannsyn
   double A = std::min(1.,  rel_prob);
  // std::cout << rel_prob << std::endl;
   //std::cout << A << std::endl;
   //aksepter eller avvis tilstanden
   double r = uniform_real(generator);
   //std::cout<< "r: " << r << std::endl;
   if (r <= A){
     //std::cout << "Flipping!" << std::endl;
     //std::cout << "Energy diff: " << energy_diff << std::endl;
     S(index_1, index_2) = -S(index_1, index_2); //flip the spin!
     make_periodic();
   }
}


//rekn ut energidifferanse
double State::delta_E(int &index_1, int &index_2)
{
  //std::cout << "Calculate energy difference" <<std::endl;
  double e_before = 0;
   e_before += S(index_1, index_2)*S(index_1+1, index_2);
   e_before += S(index_1, index_2)*S(index_1-1, index_2);
   e_before += S(index_1, index_2)*S(index_1, index_2+1);
   e_before += S(index_1, index_2)*S(index_1, index_2-1);

  return -2*e_before;
}


double State:: total_energy()
{
  //std::cout << "Kallar på midle energi" << std::endl;
  double energy = 0;
  for (int i=1;i<=L;i++){
    for (int j=1; j<= L; j++){
      energy += S(i, j)*S(i+1, j) + S(i, j)*S(i, j+1);
    }
  }
  return energy;
}


double State::total_magnetization()
{
  //std::cout << "Mag" << std::endl;
  double M = arma::accu(S.rows(1, L).cols(1, L));
 return std::abs(M);
}

double State::specific_heat_capacity()
//calculates the specific heat capacity for a state using samples from one MC cycle
{
  double c = 1/(N*T*T)*(E2_sample/N - std::pow(E_sample/N, 2));
return c;
}

double State::magnetic_susceptibility()
//calculates the magnetic susceptibility for a state using samples from one MC cycle
{
  double x = 1/(N*T)*(M2_sample/N - std::pow(M_sample/N, 2));
return x;
}
