#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include "armadillo"
#include "coshfunc.h"
#include "sampling.h"
#include "file.h"

int main(int argc, char* argv[]){
  std::string filename="../../../model_output/parameter_priors.dat";
  int lhp_samples=1000,
    ab=4,
    parameters_per_ab=6,
    parameters=ab*parameters_per_ab;

  /*Load parameter_priors information from file
    We store twice in Names[i] to rid outselves of the first field
    in the file, "UNIFORM"
   */
  arma::mat File = arma::zeros<arma::mat>(ab,2);
  std::vector<std::string> Names(parameters);
  File = load_range_file(filename,parameters,Names);

  /*Construct the generic hypercube sampling in terms of an int list, 0 to lhp_samples
    We end with a lines by lhp_s matrix. Then construct full numerical lhp sampling 
    using the input from file and the previous generic sampling matrix
  */
  arma::mat hypercube = arma::zeros<arma::mat>(lhp_samples,parameters);
  hypercube = construct_lhp_plus_chi_constraint(lhp_samples,ab,File);
  std::string printname = "../../../model_output/moments_parameters.dat";
  print_file(printname,hypercube);

  /*Construct the parameter files for each sampling 
   */
  create_parameter_prior_files(Names, hypercube);

  return 0;
}





