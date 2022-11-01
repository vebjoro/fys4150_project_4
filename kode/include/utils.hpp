#ifndef UTILS_H_
#define UTILS_H_
#include <armadillo>


//flippar ein spinn
void flip_random_spinn(arma::mat &state, int &L, double &T);

//rekn ut delta_E
double delta_E(arma::mat &state, int &index_1, int &index_2);



#endif // UTILS_H_
