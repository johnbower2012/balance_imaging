#ifndef __SAMPLING_H__
#define __SAMPLING_H__

#include <armadillo>
#include <Eigen/Dense>
#include <vector>
#include "coshfunc.h"

void ArmaToEigen(Eigen::MatrixXd &eigen, arma::mat arm);
void EigenToArma(arma::mat &arm, Eigen::MatrixXd eigen);
Eigen::MatrixXd ArmaToEigen(arma::mat &arm);
arma::mat EigenToArma(Eigen::MatrixXd &eigen);

arma::mat construct_latinhypercube_sampling(int samples, arma::mat range);
arma::mat construct_lhp_plus_chi_constraint(int samples, int ab, arma::mat range);
void LHCSampling(Eigen::MatrixXd &Sampling, int samples, int ab, Eigen::MatrixXd range);

#endif