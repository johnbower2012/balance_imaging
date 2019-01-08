#include<iostream>
#include<Eigen/Dense>
#include<string>
#include<vector>
#include<stdlib.h>
#include "system.h"
#include "analysis.h"
#include "coshfunc.h"
#include "emulator.h"
#include "mcmc.h"

int main(int argc, char* argv[]){
  //Check that usage is proper: enter the model output to use
  if(argc != 3)
    {
      printf("Usage: ./main start end\n");
      return 1;
    }

  /*************************************
   *    Write Parameters files        *
  *************************************/
  //vector to store parameter names from RangeFile
  //foldername = model BF location
  //writefilename = paramfile name, should default to parameters.dat
  //rangename = rangefile for parameters and paramnames to LCHsample for model run
  std::string 
    foldername="model_output",
    writefilename="parameters.dat",
    rangename=foldername+"/parameter_priors.dat",
    infilename=foldername+"/moments_parameters.dat",
    delimiter=" ";
  //start, finish for model runs to use
  int
    start=atoi(argv[1]),
    finish=atoi(argv[2]),
    ab=4;
  Eigen::MatrixXd Parameters;

  WriteParameterFiles(rangename, foldername, writefilename, delimiter,
		      start, finish, ab, 
		      Parameters);
  WriteFile(infilename,Parameters,delimiter);



  /*************************************
   *    Construct gab functions       *
  *************************************/
  LoadFile(infilename,Parameters,delimiter);
  WriteGABFunctions(Parameters,delimiter,ab);



  /*************************************
   *   Run Model for output           *
  *************************************/
  //Run Scott's program
  //system("cd ..;pwd; cd src/");
  


  /*************************************
   *   Load Data for analysis         *
  *************************************/
  //Model BF filenames
  std::vector<std::string> modelfilenames={"I211_J211.dat",
					   "I2212_J2212.dat",
					   "I321_J2212.dat",
					   "I321_J321.dat"};
  //EXP quark filenames
  std::vector<std::string> expfilenames={"star_pipi.dat",
					 "star_ppbar.dat",
					 "star_pK.dat",
					 "star_KK.dat"};
  //Number of model && exp files
  //data column to be selected from model/exp files
  bool
    removeValue = true;
  std::vector<Eigen::MatrixXd> 
    ModelMatrix,
    ExpMatrix;
  Eigen::VectorXd
    modeldy,
    expdy;

  LoadMEDataFiles(modelfilenames, expfilenames,
		     ModelMatrix, ExpMatrix,
		     modeldy, expdy,
		     foldername, delimiter,
		     start, finish);


  
  /*************************************
   *   Conduct model analysis         *
  *************************************/
  Eigen::MatrixXd
    ModelObs,
    ExpObs;
  //Ratio of ExpObs to Error, e.g. 0.1 --> Error=0.1*ExpObs
  double er=0.1;
  Eigen::VectorXd
    Error,
    Mean,
    EigenValues;
  Eigen::MatrixXd
    ModelTilde,
    ExpTilde,
    Covariance,
    EigenVectors,
    ModelZ,
    ExpZ;
  
  //Calculate Model & Exp Moments, calc Error
  MatrixMoments(ModelMatrix,modeldy,ModelObs);
  MatrixMoments(ExpMatrix,expdy,ExpObs);
  Error = er*ExpObs.col(0);

  //Calculate y-tilde for Model & Exp
  AverageRows(Mean,ModelObs);
  TildeFunction(ModelTilde,Mean,Error,ModelObs);
  TildeFunction(ExpTilde,Mean,Error,ExpObs);

  //Calculate Covariance of ModelTilde and related eigenvectors
  CovarianceFunction(Covariance,ModelTilde);
  EigenSolve(EigenValues,EigenVectors,Covariance);
  EigenSort(EigenValues,EigenVectors);

  //Calculate Principle Components of y-tilde
  ModelZ = ModelTilde.transpose()*EigenVectors;
  ExpZ = ExpTilde.transpose()*EigenVectors;

  //Write out PC
  WriteFile(foldername+"/expz.dat",ExpZ," ");
  WriteFile(foldername+"/modelz.dat",ModelZ," ");

  //std::cout << "Parameters:\n" << Parameters << std::endl;  
  //std::cout << "ModelObs:\n" << ModelObs << std::endl;
  //std::cout << "ExpObs:\n" << ExpObs << std::endl;
  //std::cout << "Error:\n" << Error << std::endl;
  //std::cout << "Mean:\n" << Mean << std::endl;
  //std::cout << "ModelTilde:\n" << ModelTilde << std::endl;
  //std::cout << "ExpTilde:\n" << ExpTilde << std::endl;
  //std::cout << "Covariance:\n" << Covariance << std::endl;
  std::cout << "EV:\n" << EigenValues.transpose() << std::endl;
  //std::cout << "EV:\n" << EigenVectors << std::endl;
  std::cout << "----------------\n" << std::endl;
  //std::cout << "ModelZ:\n" << ModelZ << std::endl;
  //AverageColumns(Mean,ModelZ);
  //std::cout << "ModelZMean:\n" << Mean.transpose() << std::endl;
  std::cout << "ExpZ:\n" << ExpZ << std::endl;


  /*************************************
   *   Plot ModelZ w/ Params          *
  *************************************/
  int 
    parameters = Parameters.cols(),
    observables = ModelZ.cols();
  Eigen::MatrixXd 
    plot(finish-start,parameters+observables);

  plot.block(0,0,finish-start,parameters) = Parameters;
  plot.block(0,parameters,finish-start,observables) = ModelZ;
  WriteFile("trainplot.dat",plot," ");

p
  /*************************************
   *   Conduct MCMC Analysis          *
  *************************************/
  /*
  int 
    test=finish-start;
  std::vector<std::string> 
    paramNames;
  Eigen::MatrixXd 
    Hyperparameters,
    outMatrix(test,parameters+observables),
    range,
    Beta;
  Eigen::VectorXd Width(parameters);

  //Load Hyperparameters
  LoadFile("hyperparameters.dat",Hyperparameters," ");
  //Load Range File
  LoadParamFile(rangename,paramNames,range,delimiter);
  //Create Widths
  for(int i=0;i<parameters;i++)
    {
      Width(i) = (range(i,1) - range(i,0))/50.0;
    }
  //Create & Write Beta matrix
  linearRegressionLeastSquares(ModelZ,Parameters,Beta);
  WriteFile("beta.dat",Beta," ");

  emulator emulation(Parameters, Hyperparameters, Beta);
  MCMC mcmc(ExpZ,range,Width,true);
  mcmc.setPosition();
  int Samples=10000;
  
  Eigen::MatrixXd History;
  mcmc.Run(Samples, History, emulation, ModelZ);
  Eigen::MatrixXd trace = History.block(0,0,Samples,parameters);
  WriteCSVFile("mcmctrace.csv",paramNames,trace,",");
  */

  /*************************************
   *   Conduct Fit Analysis           *
  *************************************/
  /*  
  for(int row=0;row<10;row++){
  testX = Parameters.row(row);
  Eigen::MatrixXd testY_ = ModelZ.row(row),
    PM = Parameters,
    MZ = ModelZ;
  RemoveRow(PM,row);
  RemoveRow(MZ,row);
  std::cout << "testX\n" << testX << std::endl;
  std::cout << "ModelY\n" << testY_ << std::endl;

  Eigen::MatrixXd Fit;
  emulator emulation(PM, Hyperparameters, Beta);
  outMatrix = Eigen::MatrixXd::Zero(1,parameters + observables);
  emulation.Emulate(testX, MZ, testY);
  test = 1;
  outMatrix.block(0,0,test,parameters) = testX;
  outMatrix.block(0,parameters,test,observables) = testY;
  std::cout << "WR_testY\n" << testY << std::endl;
  ComputeFit(testY_,testY,Fit);
  std::cout << "Fit:\n" << Fit << std::endl;

  emulation.Emulate_NR(testX, MZ, testY);
  test = 1;
  outMatrix.block(0,0,test,parameters) = testX;
  outMatrix.block(0,parameters,test,observables) = testY;
  std::cout << "NR_testY\n" << testY << std::endl;
  ComputeFit(testY_,testY,Fit);
  std::cout << "Fit:\n" << Fit << std::endl << std::endl;
  }
  */



  return 0;
}
