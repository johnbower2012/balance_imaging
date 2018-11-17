#ifndef __EMULATOR_CLASS_H__
#define __EMULATOR_CLASS_H__

#include<cmath>
#include<cstdlib>
#include<random>
#include<chrono>
#include "time.h"
#include <armadillo>

class emulator{
public:
  int train_points;
  int parameter_count;
  int observable_count;
  int hyperparameter_count;

  double epsilon;
  arma::vec noise;

  arma::mat X;
  arma::mat *kernel;
  arma::mat *kernel_inverse;
  arma::mat hyperparameters;
  arma::mat *H;
  arma::mat beta;

  emulator(arma::mat X, arma::mat hyperparameters, arma::mat beta, double epsilon);

  arma::mat kernel_function(arma::mat A, arma::mat B, int obs_index);
  arma::mat regression_linear_function(arma::mat X_mat, int obs_index);
  arma::mat emulate(arma::mat X_s, arma::mat Y_mat);
  arma::mat emulate_nr(arma::mat X_s, arma::vec y);
};

#endif
