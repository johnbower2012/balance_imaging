#include <armadillo>
#include "emulator_class.h"
#include "file.h"
#include "sampling.h"

int main(int argc, char* argv[]){
  int train = 1000;
  int test = 719;
  int param = 24;
  int observables=12;
  int hp = 3;
  double epsilon=1e-8;
  int INT;
  std::string fileName, betaName,hypName,dest_folder;

  arma::mat data = arma::zeros<arma::mat>(train,param+observables);
  arma::mat X = arma::zeros<arma::mat>(train,param);
  arma::mat X_s = arma::zeros<arma::mat>(test,param);
  arma::mat range = arma::zeros<arma::mat>(param,2);
  arma::mat H = arma::zeros<arma::mat>(observables,hp);
  arma::mat beta = arma::zeros<arma::mat>(observables,param+1);
  arma::mat Y = arma::zeros<arma::mat>(train,observables);

  if(argc<5){
    printf("Improper usage. Please enter 'fileName betaName hyperpName dest_folder' on same line.\n");
    exit(1);
  } else{
    fileName = argv[1];
    betaName = argv[2];
    hypName = argv[3];
    dest_folder = argv[4];
  }
  
  load_file(fileName,data);
  load_file(betaName,beta);
  load_file(hypName,H);

  for(int i=0;i<train;i++){
    for(int j=0;j<param;j++){
      X(i,j) = data(i,j);
    }
    for(int j=0;j<observables;j++){
      Y(i,j) = data(i,param+j);
    }
  }

  for(int i=0;i<param;i++){
    if(i%6==0){
      range(i,0) = 0.5;
      range(i,1) = 2.0;
    }
    else{
      range(i,0) = -2;
      range(i,1) = 2;
    }
  }

  X_s = construct_latinhypercube_sampling(test,range);

  emulator gauss(X,H,beta,epsilon);
  arma::mat output = gauss.emulate(X_s,Y);
  std::string printname = dest_folder + "test.dat";
  write_output(output,param,observables,printname);
  printname = dest_folder + "train.dat";
  write_trainset(X, Y, printname);

  return 0;
}
