#ifndef __FIT_H__
#define __FIT_H__

#include<armadillo>

arma::vec linear_regression_ls(arma::vec y, arma::mat X);
arma::mat linear_regression_ls(arma::mat Y, arma::mat X);
void add_one_column(arma::mat &X,arma::mat &X_plus_one);

#endif
