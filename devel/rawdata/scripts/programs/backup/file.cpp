#include "file.h"

void print_file(std::string outfilename, Eigen::MatrixXd matrix){
  unsigned int rows = matrix.rows();
  unsigned int cols = matrix.cols();
  std::ofstream ofile;
  ofile.open(outfilename);
  for(unsigned int i=0;i<rows;i++){
    for(unsigned int j=0;j<cols;j++){
      ofile << ' ' << matrix(i,j);
    }
    ofile << '\n';
  }
  ofile.close();
}
void load_file(std::string infilename, Eigen::MatrixXd &matrix, unsigned int parameters, unsigned int samples){
  matrix = Eigen::MatrixXd::Zero(samples,parameters);
  std::ifstream ifile;
  ifile.open(infilename);
  for(unsigned int i=0;i<samples;i++){
    for(unsigned int j=0;j<parameters;j++){
      ifile >> matrix(i,j);
    }
  }
  ifile.close();
}

/////////////////////////////

void load_file(int files, int lines, int runs, std::string *&infilename, arma::vec &delY_vec, arma::mat *&val_matrix){
  std::ifstream ifile;
  if(val_matrix!=nullptr){
    delete[] val_matrix;
  }
  val_matrix = new arma::mat[files];
  int i, j, k;

  for(i=0;i<files;i++){
    val_matrix[i] = arma::zeros<arma::mat>(runs,lines);
    ifile.open(infilename[i]);
    printf("+++++ %s LOADED +++++\n", infilename[i].c_str());
    for(j=0;j<lines;j++){
      for(k=0;k<runs+1;k++){
	if(k==0){
	  ifile >> delY_vec(j);
	} else{
	  ifile >> val_matrix[i](k-1,j);
	}
      }
    }
    ifile.close();
  }
}
void print_file(std::string outfilename, std::string title, arma::vec vector){
  int elem = vector.n_elem;
  std::ofstream ofile;
  ofile.open(outfilename);
  ofile << title.c_str() << '\n';
  for(int i=0;i<elem;i++){
    ofile << ' ' << vector(i);
  }
  ofile << '\n';
  ofile.close();
}
void print_file(std::string outfilename, arma::vec vector){
  int elem = vector.n_elem;
  std::ofstream ofile;
  ofile.open(outfilename);
  for(int i=0;i<elem;i++){
    ofile << ' ' << vector(i);
  }
  ofile << '\n';
  ofile.close();
}
void print_file(std::string outfilename, std::string title, arma::mat matrix){
  int rows = matrix.n_rows;
  int cols = matrix.n_cols;
  std::ofstream ofile;
  ofile.open(outfilename);
  ofile << title.c_str() << '\n';
  for(int i=0;i<rows;i++){
    for(int j=0;j<cols;j++){
      ofile << ' ' << matrix(i,j);
    }
    ofile << '\n';
  }
  ofile.close();
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


void create_parameter_prior_files(std::vector<std::string> Names, arma::mat hypercube){
  std::ofstream ofile;
  int lhp_samples = hypercube.n_rows;
  int parameters = hypercube.n_cols;
  for(int i=0;i<lhp_samples;i++){
    char fn[50];
    sprintf(fn,"../../../model_output/run%04d/parameters.dat",i);
    ofile.open(fn);
    for(int j=0;j<parameters;j++){
      ofile << Names[j] << " " << hypercube(i,j) << '\n';
    }
    ofile.close();
  }
}
