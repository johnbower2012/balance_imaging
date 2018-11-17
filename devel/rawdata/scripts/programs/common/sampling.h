#ifndef __SAMPLING_H__
#define __SAMPLING_H__

#include <armadillo>
#include <Eigen/Dense>
#include "coshfunc.h"

arma::mat construct_latinhypercube_sampling(int samples, arma::mat range);
arma::mat construct_lhp_plus_chi_constraint(int samples, int ab, arma::mat range);

#endif
