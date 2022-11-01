#include "utils.hpp"
#include <armadillo>
#include <random>
#include <chrono>

//flippar ein spinn
void flip_random_spinn(arma::mat &state, int &L, double &T)
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
  state(index_1, index_2) = -state(index_1, index_2);
   double energy_diff = delta_E(state, index_1, index_2);
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
     state(index_1, index_2) = -state(index_1, index_2); //flip it back
   }

}




//rekn ut energidifferanse
double delta_E(arma::mat &state, int &index_1, int &index_2)
{
  //std::cout << "Calculate energy difference" <<std::endl;
  double D_inter_E = 0;
  D_inter_E += state(index_1, index_2)*state(index_1+1, index_2);
  D_inter_E += state(index_1, index_2)*state(index_1-1, index_2);
   D_inter_E += state(index_1, index_2)*state(index_1, index_2+1);
   D_inter_E += state(index_1, index_2)*state(index_1, index_2-1);
   D_inter_E = 2*D_inter_E;

  return D_inter_E;
}
