#ifndef __EMULATOR_CLASS_H__
#define __EMULATOR_CLASS_H__

#include<cmath>
#include<cstdlib>
#include<random>
#include<chrono>
#include<fstream>
#include "time.h"
#include "armadillo"

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

arma::mat construct_latinhypercube_sampling(int samples, arma::mat range);
void write_output(arma::mat output_mat, int param, int obs, std::string outfilename);
void write_trainset(arma::mat X_mat, arma::mat Y_mat, std::string outfilename);
void load_data_file(std::string fileName, arma::mat &X, arma::mat &Y);
void load_beta_file(std::string betaName, arma::mat &beta);
void load_file(std::string fileName, arma::mat &file);
void write_file(std::string fileName, arma::mat &file);

#endif