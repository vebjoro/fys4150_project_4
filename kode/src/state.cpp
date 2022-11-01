#include "state.hpp"
#include <armadillo>
#include <random>
#include <chrono>



//konstruktør
State::State(int size, double temp){
  L = size;
  T = temp;

  int N = L*L;
  arma::mat S = arma::mat(L+2, L+2);
  S = arma::sign(S.randu() - 0.5);
  S.replace(0, 1);

  //periodisitet
  arma::vec p_right = S.col(1);
  arma::vec p_left = S.col(L);
  S.col(0) = p_left;
  S.col(L+1) = p_right;

   arma::rowvec p_top = S.row(1);
   arma::rowvec p_bottom = S.row(L);
   S.row(0) = p_bottom;
  S.row(L+1) = p_top;


  S(0, 0) = 0;
  S(L+1, L+1) = 0;
  S(L+1, 0) = 0;
  S(0, L+1) = 0;
  std::cout << S << std::endl;
}



//flippar ein spinn
void State::flip_random_spinn()
{
  //std::cout << "Calling flip spinn" << std::endl;

  unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
  //Construct Mersenne Twister and seed it
  std::mt19937 generator;
  std::uniform_int_distribution<int> uniform_dist(1, L);
  std::uniform_real_distribution<> uniform_real(0, 1);
  generator.seed(seed);

  //brukar generatoren til å velge ein tilfeldig spinn
  int index_1  = uniform_dist(generator); // korleis kalle på arma-element med ein index?
  int index_2  = uniform_dist(generator);
  //std::cout << "Indeksar: " << index_1  << ", " << index_2 << std::endl;

  //flip the spin!
  S(index_1, index_2) = -S(index_1, index_2);
   double energy_diff = delta_E(index_1, index_2);
   //std::cout << energy_diff << std::endl;
   double rel_prob = std::exp(-energy_diff/T); //double check this!


   //godkjenningssannsyn
   double A = std::min(1.,  rel_prob);



   //aksepter eller avvis tilstanden
   double r = uniform_real(generator);
   //std::cout<< "r: " << r << std::endl;
   if (r <= A){
     //std::cout << "Flipping!" << std::endl; do nothing, spin is already flipped
   }
   else{
     S(index_1, index_2) = -S(index_1, index_2); //flip it back
   }

}




//rekn ut energidifferanse
double State::delta_E(int &index_1, int &index_2)
{
  //std::cout << "Calculate energy difference" <<std::endl;
  double D_inter_E = 0;
  D_inter_E += S(index_1, index_2)*S(index_1+1, index_2);
  D_inter_E += S(index_1, index_2)*S(index_1-1, index_2);
   D_inter_E += S(index_1, index_2)*S(index_1, index_2+1);
   D_inter_E += S(index_1, index_2)*S(index_1, index_2-1);
   D_inter_E = 2*D_inter_E;

  return D_inter_E;
}
