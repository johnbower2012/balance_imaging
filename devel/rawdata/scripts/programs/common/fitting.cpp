#include "fitting.h"

arma::mat linear_regression_ls(arma::mat Y, arma::mat X){
  int points = Y.n_rows;
  int parameters = X.n_cols-1;
  int observables = Y.n_cols;
  arma::mat temp = arma::zeros<arma::mat>(points,2);

  arma::mat Beta = arma::zeros<arma::mat>(observables,1+parameters);
  arma::vec beta = arma::zeros<arma::vec>(1+parameters);
  arma::vec y = arma::zeros<arma::vec>(points);

  temp = X.t()*X;
  temp = temp.i();
  for(int j=0;j<observables;j++){
    y = Y.col(j);
    beta = temp*X.t()*y;
    Beta.row(j) = beta.t();
  }

  return Beta;
}
void add_one_column(arma::mat &X,arma::mat &X_plus_one){
  int rows = X.n_rows,
    cols = X.n_cols;
  X_plus_one = arma::zeros<arma::mat>(rows,cols+1);
  for(int i=0;i<rows;i++){
    X_plus_one(i,0) = 1.0;
    for(int j=0;j<cols;j++){
      X_plus_one(i,j+1) = X(i,j);
    }
  }
}
