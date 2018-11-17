#ifndef __OBSERVABLES_CLASS_H__
#define __OBSERVABLES_CLASS_H__

#include<fstream>
#include<cmath>
#include "armadillo"

arma::vec average_columns(arma::mat input);
void print_fractional_sum(arma::vec vector);

void tilde_function(const arma::mat &matrix, const arma::vec &error, arma::vec &mean, arma::mat &tilde);
void covariance_function(const arma::mat &matrix, arma::mat &covariance);
void sort_eigen_function(arma::vec &eigval, arma::mat &eigvec);

double zeroth_moment(const arma::mat &function);
double first_moment(const arma::mat &function);
double second_moment(const arma::mat &function);
void obs_matrix_moments(int files, int obs_file, arma::mat *&val_matrix, const arma::vec &delY_vec, arma::mat &obs_matrix);


#endif
