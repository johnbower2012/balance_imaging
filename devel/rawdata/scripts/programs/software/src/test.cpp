#include<iostream>
#include<Eigen/Dense>
#include<string>
#include<vector>
#include<stdlib.h>
#include "system.h"
#include "analysis.h"

int main(int argc, char* argv[]){

  std::vector<std::string> modelfilenames={"I211_J211.dat","I2212_J2212.dat","I321_J2212.dat","I321_J321.dat"},
    expfilenames={"star_pipi.dat","star_ppbar.dat","star_pK.dat","star_KK.dat"};

  std::string foldername="model_output",
    writefilename="parameters.dat",
    rangename=foldername+"/parameter_priors.dat",
    delimiter=" ";

  int start=atoi(argv[1]),
    finish=atoi(argv[2]),
    column=1,
    ab=4,
    modelfiles = modelfilenames.size(),
    expfiles = expfilenames.size();

  std::vector<Eigen::MatrixXd> ModelMatrix(modelfiles),
    ExpMatrix(expfiles);
  Eigen::MatrixXd Dy,ModelObs,ExpObs,Parameters;
  Eigen::VectorXd modeldy, expdy;
  LoadDataFile(foldername, modelfilenames[0], delimiter, 0, 1, 0, Dy);
  modeldy = Dy.col(0);
  LoadDataFile(foldername, expfilenames[0], delimiter, 0, 1, 0, Dy);
  expdy = Dy.col(0);

  //WriteParameterFiles(rangename, foldername, writefilename, delimiter, start, finish, ab,Parameters);
  //system("cd ..;pwd; cd src/");
  LoadDataFiles(foldername, modelfilenames, delimiter, start, finish, column, ModelMatrix);
  LoadDataFiles(foldername, expfilenames, delimiter, 0, 1, column, ExpMatrix);
  MatrixMoments(ModelMatrix,modeldy,ModelObs);
  MatrixMoments(ExpMatrix,expdy,ExpObs);
  
  double er=0.1;
  Eigen::VectorXd Error = er*ExpObs.col(0),
    Mean,
    EigenValues;
  Eigen::MatrixXd ModelTilde,
    ExpTilde,
    Covariance,
    EigenVectors,
    ModelZ,
    ExpZ;

  AverageRows(Mean,ModelObs);
  TildeFunction(ModelTilde,Mean,Error,ModelObs);
  TildeFunction(ExpTilde,Mean,Error,ExpObs);

  CovarianceFunction(Covariance,ModelObs);
  EigenSolve(EigenValues,EigenVectors,Covariance);
  EigenSort(EigenValues,EigenVectors);
  ModelZ = ModelObs.transpose()*EigenVectors;
  ExpZ = ExpObs.transpose()*EigenVectors;

  //std::cout << "Parameters:\n" << Parameters << std::endl;  
  //std::cout << "ModelObs:\n" << ModelObs << std::endl;
  //std::cout << "ExpObs:\n" << ExpObs << std::endl;
  //std::cout << "Error:\n" << Error << std::endl;
  //std::cout << "Mean:\n" << Mean << std::endl;
  //std::cout << "ModelTilde:\n" << ModelTilde << std::endl;
  std::cout << "ExpTilde:\n" << ExpTilde << std::endl;
  std::cout << "Covariance:\n" << Covariance << std::endl;
  std::cout << "EV:\n" << EigenValues << std::endl;
  std::cout << "EV:\n" << EigenVectors << std::endl;
  std::cout << "----------------\n" << std::endl;
  //std::cout << "ModelZ:\n" << ModelZ << std::endl;
  AverageColumns(Mean,ModelZ);
  std::cout << "Mean:\n" << Mean.transpose() << std::endl;
  std::cout << "ExpZ:\n" << ExpZ << std::endl;
  

  return 0;
}