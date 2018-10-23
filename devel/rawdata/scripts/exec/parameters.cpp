#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include "armadillo"
#include "../lib/coshfunc.cpp"

void print_file(std::string outfilename, arma::mat matrix);
arma::mat load_range_file(std::string filename, int parameters, std::vector<std::string> &Names);
arma::mat construct_latinhypercube_sampling(int samples, arma::mat range);
arma::mat construct_lhp_plus_chi_constraint(int samples, int ab, arma::mat range);
arma::vec FCoshIntegral(int nmax);
void create_parameter_prior_files(std::vector<std::string> Names, arma::mat hypercube);

int main(int argc, char* argv[]){
  std::string filename;
  int lines,lhp_samples,ab=4;
  if(argc!=4){
    std::cout << "Usage: Enter also 'fn lines lhp_s' the same line.\n";
    exit(1);
  }
  else{
    filename=argv[1];
    lines=atoi(argv[2]);
    lhp_samples=atoi(argv[3]);
  }

  /*Load parameter_priors information from file
    We store twice in Names[i] to rid outselves of the first field
    in the file, "UNIFORM"
   */
  arma::mat File = arma::zeros<arma::mat>(lines,2);
  std::vector<std::string> Names(lines);
  File = load_range_file(filename,lines,Names);

  /*Construct the generic hypercube sampling in terms of an int list, 0 to lhp_samples
    We end with a lines by lhp_s matrix. Then construct full numerical lhp sampling 
    using the input from file and the previous generic sampling matrix
  */
  arma::mat hypercube = arma::zeros<arma::mat>(lhp_samples,lines);
  hypercube = construct_lhp_plus_chi_constraint(lhp_samples,ab,File);
  std::string printname = "../model_output/moments_parameters.dat";
  print_file(printname,hypercube);

  /*Construct the parameter files for each sampling 
   */
  create_parameter_prior_files(Names, hypercube);

  return 0;
}


void print_file(std::string outfilename, arma::mat matrix){
  int rows = matrix.n_rows;
  int cols = matrix.n_cols;
  std::ofstream ofile;
  ofile.open(outfilename);
  for(int i=0;i<rows;i++){
    for(int j=0;j<cols;j++){
      ofile << ' ' << matrix(i,j);
    }
    ofile << '\n';
  }
  ofile.close();
}
arma::mat load_range_file(std::string filename, int parameters, std::vector<std::string> &Names){
  std::ifstream ifile;
  arma::mat File = arma::zeros<arma::mat>(parameters,2);
  Names.resize(parameters);
  ifile.open(filename);
  for(int i=0;i<parameters;i++){
    ifile >> Names[i]; ifile >> Names[i];
    ifile >> File(i,0);
    ifile >> File(i,1);
  }
  ifile.close();
  return File;
}
arma::mat construct_latinhypercube_sampling(int samples, arma::mat range){
  int parameters = range.n_rows;
  arma::mat hypercube = arma::zeros<arma::mat>(samples,parameters);
  arma::vec hyperlist = arma::linspace<arma::vec>(0,samples-1,samples);
  for(int i=0;i<parameters;i++){
    hyperlist=shuffle(hyperlist);
    hypercube.col(i) = hyperlist;
  }

  float init,final,dx;
  for(int i=0;i<parameters;i++){
    init = range(i,0);
    final = range(i,1);
    dx = (final-init)/(samples-1);
    hyperlist = hypercube.col(i);
    hyperlist = init + dx*hyperlist;
    hypercube.col(i) = hyperlist;
  }
  return hypercube;
}  
arma::mat construct_lhp_plus_chi_constraint(int samples, int ab, arma::mat range){
  int parameters = range.n_rows;
  int nmax = parameters/ab - 2;
  CCosh dist(nmax);
  Eigen::VectorXd G = Eigen::VectorXd::Zero(ab);

  arma::mat hypercube = arma::zeros<arma::mat>(samples,parameters);
  arma::vec hyperlist = arma::linspace<arma::vec>(0,samples-1,samples);
  arma::mat chi = {{ 7.2728e+02, -2.2654e+02, -8.3627e+01},
		   {-2.2654e+02,  7.2395e+02, -8.2073e+01},
		   {-8.3627e+01, -8.2073e+01,  3.0778e+02}};
  arma::vec Chi(4);
  Chi(0) = chi(0,0);
  Chi(1) = chi(0,1);
  Chi(2) = chi(0,2);
  Chi(3) = chi(2,2);
  for(int i=0;i<parameters;i++){
    hyperlist=shuffle(hyperlist);
    hypercube.col(i) = hyperlist;
  }

  float init,final,dx;
  for(int i=0;i<parameters;i++){
    init = range(i,0);
    final = range(i,1);
    dx = (final-init)/(samples-1);
    hyperlist = hypercube.col(i);
    hyperlist = init + dx*hyperlist;
    hypercube.col(i) = hyperlist;
  }

  for(int i=0;i<samples;i++){
    for(int j=0;j<ab;j++){
      hypercube(i,j*(nmax+2)+1) = 1.0/hypercube(i,j*(nmax+2));
      for(int k=1;k<nmax+1;k++){
	hypercube(i,j*(nmax+2)+1) -= hypercube(i,j*(nmax+2)+k+1)*dist.Z(k);
      }
      hypercube(i,j*(nmax+2)+1) /= (dist.Z(0));
    }
  }

  /*
  for(int i=0;i<samples;i++){
    G = Eigen::VectorXd::Zero(ab);
    for(int j=0;j<ab;j++){
      for(int k=0;k<nmax+1;k++){
	G(j) += hypercube(i,j*(nmax+2)+k+1)*dist.Z(k);
      }
      for(int k=0;k<nmax+1;k++){
	hypercube(i,j*(nmax+2)+k+1) = hypercube(i,j*(nmax+2)+k+1)*Chi(j)/G(j);
      }
    }
  }
  */
  return hypercube;
}  

void create_parameter_prior_files(std::vector<std::string> Names, arma::mat hypercube){
  std::ofstream ofile;
  int lhp_samples = hypercube.n_rows;
  int parameters = hypercube.n_cols;
  for(int i=0;i<lhp_samples;i++){
    char fn[50];
    sprintf(fn,"../model_output/run%04d/parameters.dat",i);
    ofile.open(fn);
    for(int j=0;j<parameters;j++){
      ofile << Names[j] << " " << hypercube(i,j) << '\n';
    }
    ofile.close();
  }
}







