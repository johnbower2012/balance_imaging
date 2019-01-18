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

int main(void)
{
  std::vector<std::string> 
    Name;
  std::string 
    foldername="model_output",
    filename="posterior.dat",
    writefilename="parameters.dat",
    rangename=foldername+"/parameter_priors.dat",
    infilename=foldername+"/moments_parameters.dat",
    delimiter=" ";
  //start, finish for model runs to use
  int
    start=0,
    finish=20,
    ab=4;
  Eigen::MatrixXd 
    matrix,
    range,
    mean,
    std;
  LoadParamFile(rangename,Name,range,delimiter);
  WritePosteriorParameterFiles(foldername,filename,Name,writefilename,delimiter,start,finish);
  std::string cmd="for((i=0;i<20;i++)); do fn=$(printf 'model_output/run%04d/' $i); cp -v model_output/fixed_parameters.dat ${fn}; done";
  system(cmd.c_str());

}