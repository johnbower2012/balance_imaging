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
  std::string 
    foldername="model_output",
    writefilename="parameters.dat",
    rangename=foldername+"/parameter_priors.dat",
    infilename=foldername+"/moments_parameters.dat",
    delimiter=" ";
  //start, finish for model runs to use
  int
    start=0,
    finish=1000,
    ab=4;
  Eigen::MatrixXd Parameters;

  WriteParameterFiles(rangename, foldername, writefilename, delimiter,
		      start, finish, ab, 
		      Parameters);
  WriteFile(infilename,Parameters,delimiter);
  std::string cmd = "for((i=0;i<1000;i++)); do fn=$(printf "run%04d/" $i); cp -v fixed_parameters.dat $fn; done";
  system(cmd);
}
