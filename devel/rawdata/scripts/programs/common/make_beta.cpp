#include<iostream>
#include<armadillo>
#include<string>
#include "fitting.h"
#include "file.h"

int main(int argc, char* argv[]){
  int points=1000,dimensions=4,observables=12;

  arma::mat x = arma::zeros<arma::mat>(points,dimensions);
  arma::mat X = arma::zeros<arma::mat>(points,dimensions+1);
  arma::mat Y = arma::zeros<arma::mat>(points,observables);
  arma::vec y = arma::zeros<arma::vec>(points);
  arma::vec result = arma::zeros<arma::vec>(points);
  arma::vec beta = arma::zeros<arma::vec>(points);
  std::string plotname;
  std::string dest_folder;

  if(argc<3){
    printf("Improper usage. Enter also 'plotname dest_folder' on same line.\n");
    exit(1);
  } else{
    plotname = argv[1];
    dest_folder = argv[2];
  }
  
  load_data_file(plotname,x,Y);
  add_one_column(x,X);

  arma::mat Beta = linear_regression_ls(Y,X);

  std::string destfile = dest_folder + "beta.dat";
  write_file(destfile,Beta);

  return 0;
}


