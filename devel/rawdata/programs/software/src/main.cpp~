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
#include "parametermap.h"

int main(int argc, char* argv[]){
  CParameterMap info;
  string runname = "run.dat";
  info.ReadParsFromFile(runname);
  info.PrintPars();
  int
    start = info.getI("START",0),
    finish = info.getI("FINISH",200),
    action = info.getI("ACTION",1),
    ab = info.getI("QUARK_PAIRS",4);

  /*************************************
   *    Write Parameters files        *
  *************************************/
  std::string 
    infoldername=info.getS("MODEL_FOLDER","model_folder"),          //infoldername = model BF location
    outfoldername=info.getS("OUTPUT_FOLDER","stat_output"),         //outfoldername = write destination
    writefilename=info.getS("PARAM_FILE","parameters.dat"),         //writefilename = paramfile name, should default to parameters.dat
    rangename=infoldername+"/"+info.getS("PARAM_RANGE_FILE","parameter_priors.dat"),    //rangename = rangefile for parameters and paramnames to LCHsample for model run
    infilename=outfoldername+"/"+info.getS("PARAMS_FILE","parameters_lch.dat"), //output parameters used for model runs
    delimiter=" ",   //delimiter for writing and reading files
    gabname=outfoldername+"/"+info.getS("QUARK_FUNC_FILE","gabfunctions.dat");  //write prior gab functions, in order eta, uu0, ud0, us0, ss0, uu1, ud1, ...

  Eigen::MatrixXd Parameters;

  WriteParameterFiles(rangename, infoldername, writefilename, delimiter,
		      start, finish, ab, 
		      Parameters);

  WriteFile(infilename,Parameters,delimiter);

  /*************************************
   *    Construct gab functions       *
  *************************************/
  WriteGABFunctions(gabname,Parameters,delimiter,ab);

  /*************************************
   *   Run Model for output           *
  *************************************/
  //Run Scott's program
  //std::string cmd = "./run.sh "+std::to_string(start)+" "+std::to_string(finish);
  //std::cout << cmd << std::endl;
  //system(cmd.c_string());
  

  /*************************************
   *   Conduct model analysis         *
  *************************************/
  Eigen::MatrixXd 
    ModelZ,
    ExpZ;
  ConductModelAnalysis(infoldername, outfoldername, delimiter, start, finish, ModelZ, ExpZ);

  /*************************************
   *   Plot ModelZ w/ Params          *
  *************************************/
  //for(int i=ab-1;i>-1;i--)
  //  {
  //    RemoveColumn(Parameters,1+i*3);
  //  }

  int 
    parameters=Parameters.cols(),
    observables=ModelZ.cols(),
    samples=finish-start;

  Eigen::MatrixXd 
    plot(finish-start,parameters+observables);
  Eigen::MatrixXd MinMax, ScaledParam;

  ScaleMatrixColumnsUniform(Parameters,MinMax,ScaledParam);

  plot.block(0,parameters,finish-start,observables) = ModelZ;

  plot.block(0,0,finish-start,parameters) = Parameters;

  WriteFile(outfoldername+"/trainplot.dat",plot," ");
  plot.block(0,0,finish-start,parameters) = ScaledParam;
  WriteFile(outfoldername+"/scaledtrainplot.dat",plot," ");
  WriteFile(outfoldername+"/scaledmoments_parameters.dat",ScaledParam," ");

  /*************************************
   *   Conduct MCMC Analysis          *
  *************************************/
  if(action==1)
    {
      int 
	test=finish-start;
      double
	fraction=1.0/300.0;
      std::vector<std::string> 
	paramNames;
  Eigen::MatrixXd 
    Hyperparameters,
    outMatrix(test,parameters+observables),
    range,
    Beta;
  Eigen::VectorXd
    Width;

  //Load Hyperparameters
  //LoadFile("hyperparameters.dat",Hyperparameters," ");
  ConstructHyperparameters(ModelZ,Hyperparameters);
  WriteFile(outfoldername+"/hyperparameters.dat",Hyperparameters,delimiter);
  //Load Range File
  LoadParamFile(rangename,paramNames,range,delimiter);

  for(int i=ab-1;i>-1;i--)
    {
      RemoveRow(MinMax,i*3+1);
      RemoveRow(range,i*3+1);
      RemoveColumn(ScaledParam,i*3+1);
    }

  //Create Widths
  ConstructWidths(Width,range,parameters,fraction);
  
  linearRegressionLeastSquares(ModelZ,ScaledParam,Beta);
  WriteFile(outfoldername+"/scaledbeta.dat",Beta," ");

  int row=50;
  ExpZ = ModelZ.row(row);
  std::cout << row << " row: " << ScaledParam.row(row) << " " << ExpZ << std::endl;
  
  Eigen::MatrixXd Goal = Eigen::MatrixXd::Zero(1,parameters+observables);
  Goal.block(0,0,1,parameters) = Parameters.row(row);
  Goal.block(0,parameters,1,observables) = ExpZ;
  WriteFile(outfoldername+"/target.dat",Goal,delimiter);

  parameters = ScaledParam.cols();
  Goal = Eigen::MatrixXd::Zero(1,parameters+observables);
  Goal.block(0,0,1,parameters) = ScaledParam.row(row);
  WriteFile(outfoldername+"/scaledtarget.dat",Goal,delimiter);

  emulator emulation(ScaledParam, Hyperparameters, Beta);
  bool regression=false;
  MCMC mcmc(ExpZ,range,Width,regression);
  mcmc.setPosition();
  int NSamples=100000;

  Eigen::MatrixXd History;
  printf("Step width fraction: %f\n",fraction);
  mcmc.Run(NSamples, History, emulation, ModelZ, MinMax);
  Eigen::MatrixXd 
    trace = History.block(0,0,NSamples,parameters),
    ntrace,
    posterior;

  std::string posteriorname = outfoldername+"/scaledposterior.dat";
  WriteFile(outfoldername+"/minmax.dat",MinMax,delimiter);

  Extract5(ntrace,trace);
  WriteCSVFile(outfoldername+"/scaledmcmctrace.csv",paramNames,ntrace,",");

  ExtractOnly20(posterior,ntrace);
  printf("Writing posterior.dat files...\n");  
  WriteFile(posteriorname,posterior,delimiter);
  for(int i=0,count=posterior.rows();i<count;i++)
    {
      for(int j=0;j<parameters;j++)
	{
	  posterior(i,j) = (posterior(i,j)*(MinMax(j,1) - MinMax(j,0)) + MinMax(j,0));
	}
    }
  for(int i=0,count=ntrace.rows();i<count;i++)
    {
      for(int j=0;j<parameters;j++)
	{
	  ntrace(i,j) = (ntrace(i,j)*(MinMax(j,1) - MinMax(j,0)) + MinMax(j,0));
	}
    }
  
  parameters += 4;
  Eigen::MatrixXd
    fullposterior = Eigen::MatrixXd::Zero(posterior.rows(),parameters),
    fullntrace = Eigen::MatrixXd::Zero(ntrace.rows(),parameters);
  int nmax=1;

  CCosh dist(5);
  for(int i=0;i<posterior.rows();i++){
    for(int j=0;j<ab;j++){
      fullposterior(i,j*(nmax+2)) = posterior(i,j*(nmax+1));
      fullposterior(i,j*(nmax+2)+2) = posterior(i,j*(nmax+1)+1);
      fullposterior(i,j*(nmax+2)+1) = 1.0/fullposterior(i,j*(nmax+2));
      for(int k=1;k<nmax+1;k++){
	fullposterior(i,j*(nmax+2)+1) -= fullposterior(i,j*(nmax+2)+k+1)*dist.Z(k);
      }
      fullposterior(i,j*(nmax+2)+1) /= (dist.Z(0));
    }
  }      

  posteriorname = outfoldername + "/posterior.dat";
  ExtractOnly20(posterior,ntrace);
  WriteFile(posteriorname,posterior,delimiter);

  gabname = outfoldername + "/posteriorgabfunctions.dat";
  WriteGABFunctions(gabname,fullposterior,delimiter,ab);

  WriteCSVFile(outfoldername+"/mcmctrace.csv",paramNames,ntrace,",");
  
  Eigen::VectorXd
    avg = Eigen::VectorXd::Zero(parameters-4);
  AverageColumns(avg,ntrace);
  Eigen::MatrixXd
    Avg = Eigen::MatrixXd::Zero(1,parameters-4);
  Avg.row(0) = avg;
  WriteFile(outfoldername+"/mcmcmean.dat",Avg,delimiter);
}


  /*************************************
   *   Conduct Fit Analysis           *
  *************************************/
  Eigen::MatrixXd
    print = Eigen::MatrixXd::Zero(finish-start,parameters+observables);

  if(action==2)
    {
      for(int i=ab-1;i>-1;i--)
	{
	  RemoveColumn(ScaledParam,1+i*3);
	}
      parameters = ScaledParam.cols();
      for(int i=0;i<finish-start;i++){
	int 
	  test=finish-start;
	std::vector<std::string> 
	  paramNames;
	Eigen::MatrixXd 
	  Hyperparameters,
	  outMatrix(test,parameters+observables),
	  Beta;
	Eigen::VectorXd
	  Width;
  
      //Load Hyperparameters
      //LoadFile("hyperparameters.dat",Hyperparameters," ");
	int row=i;
	samples = ModelZ.rows();
      Eigen::MatrixXd
	testX = ScaledParam.row(row),
	Fit,

	testY = ModelZ.row(row),
	trainY = ModelZ,
	trainX= ScaledParam;

      RemoveRow(trainX,row);
      RemoveRow(trainY,row);

      ConstructHyperparameters(trainY,Hyperparameters);  
      linearRegressionLeastSquares(trainY,trainX,Beta);

      emulator emulation(trainX, Hyperparameters, Beta);
      emulation.Emulate_NR(testX,trainY,outMatrix);

      print.block(row,0,1,parameters) = testX;
      print.block(row,parameters,1,observables) = outMatrix;
      std::cout << i << ":" << std::endl;
      std::cout << testY << std::endl;
      std::cout << outMatrix << std::endl;

      }  
  WriteFile("outMatrix.dat",print," ");
    }


  return 0;
}
